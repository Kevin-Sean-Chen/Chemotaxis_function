% mGLM
%%% mixture GLM infernece with ground truth generative model
%% define variables
lt = 50000; % length of simulation / data
nB = 4;  % number of basis function for the kernel
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3); % basis function
%%% true params
beta = 1.; % nonlinear parameter
alpha_h = [-4:-1]*.1;  % angle history kernel coefficient
alpha_dc = [1,4,-2,-1]*1;  % dC kernel coefficient
alpha_dcp = [-4:-1]*.01;  % dCp kernel coefficient
base = 0;  %baseline
kappa_turn = 5;  % turning angle variance
kappa_wv = 10;  % weather-vaning angle variance
gamma = 0.3;  % turning mixture parameter (weight on uniform)
A = 0.5;  % maximum turn probability

%% generate data
K_h = fliplr(alpha_h*cosBasis');  % dth kernel
K_dc = fliplr(alpha_dc*cosBasis');  % dC kernel
K_dcp = fliplr(alpha_dcp*cosBasis');  % dCp kernel
lc = 10;  %length of smoothing
dC = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dC stimulus vector
dCp = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dCp stimulus vector
dth = zeros(1,lt);
turns = zeros(1,lt);
F = dth*0;
pad = length(K_h);
for tt=pad+1:lt
    F(tt) = dC(tt-pad+1:tt)*K_dc' + abs(dth(tt-pad:tt-1))*K_h';  % linear filtering
    turns(tt) = choice(NL(F(tt)+base,A));  % nonlinearity and binary choice
    if rand<gamma
        mix_th = circ_vmrnd(0,0.,1)-pi;
    else
        mix_th = circ_vmrnd(pi,kappa_turn,1);
    end
    dth(tt) = turns(tt)*mix_th + (1-turns(tt))*circ_vmrnd(dCp(tt-pad+1:tt)*K_dcp',kappa_wv,1);  % angle drawn from mixture of von Mesis
          % wrapToPi((1)*circ_vmrnd(pi,kappa_turn,1)+gamma*(circ_vmrnd(0,0.,1)-pi))
end

figure; hist(dth,100); xlim([-pi,pi])
%% MLE inference
lfun = @(x)nLL(x, dth, dCp, dC, cosBasis, 0.5);  % objective function
opts = optimset('display','iter');
LB = [ones(1,12)*-10, 0, 0, 0, 0]*1;
UB = [ones(1,12)*10, 20, 20, 1, 1]*1;
prs0 = [alpha_h, alpha_dc, alpha_dcp, kappa_turn, kappa_wv, gamma, A];%
num_par = length(prs0);
prs0 = prs0 + rand(1,num_par)*0.2.*prs0;
[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);  % constrained optimization
% [x,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(lfun, prs0, opts);
fval

%% evaluation
beta_rec = x(1);
base_rec = x(end);
K_h_rec = fliplr(x(1:4)*cosBasis');
K_dc_rec = fliplr(x(5:8)*cosBasis');
K_dcp_rec = fliplr(x(9:12)*cosBasis');
figure()
subplot(131)
plot(K_h_rec); hold on; plot(K_h,'--');
subplot(132)
plot(K_dc_rec); hold on; plot(K_dc,'--');
subplot(133)
plot(K_dcp_rec); hold on; plot(K_dcp,'--');

%% evaluate density
K1 = x(13); K2 = x(14); gamma = x(15); A_ = x(16);
filt_ddc = conv_kernel(dC, fliplr(K_dc_rec));
filt_dth = conv_kernel(abs(dth), fliplr(K_h_rec));
dc_dth = filt_ddc + 1*filt_dth;
Pturns = NL(dc_dth,A_); %1./ (1 + exp( -(dc_dth) )) + 0;
n_brw = sum(Pturns)*1;
n_wv = sum(1-Pturns);
p_z = n_brw + n_wv;
p_brw = n_brw/p_z;
p_wv = n_wv/p_z;
filt_dcp = conv_kernel(dCp, fliplr(K_dcp_rec));
figure;
nbins = 100;
hh = histogram(dth, nbins, 'Normalization', 'probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
% [aa,bb] = hist((dth - filt_dcp)*1 , 1000);
bb = hh.BinEdges(1:end-1);
scal = sum(1/(2*pi*besseli(0,K2)) * exp(K2*cos( bb )) * p_wv + (1/(2*pi*besseli(0,K1)) * exp(K1*cos( bb-pi ))*(1-gamma) + (gamma)/(2*pi)) * p_brw); 
plot( bb, 1/(2*pi*besseli(0,K2)) * exp(K2*cos( bb )) * p_wv / scal , 'b'); hold on
plot( bb, ((1/(2*pi*besseli(0,K1)) * exp(K1*cos( bb-pi ))*(1-gamma) + (gamma)/(2*pi)) * p_brw)/scal ,'r'); xlim([-pi,pi])
title('von Mises for \delta C^{\perp}')

%% evaluate decision function
figure
xx = linspace(-5,5,50);
kn = norm(K_dc);
kn_ = norm(K_dc_rec);
plot(xx, A./(1+exp(-kn*xx))); hold on
plot(xx, A_./(1+exp(-kn_*xx)));

%% functions
%%% nonlinear function
function [P] = NL(F,A)
    beta = 1;
    P = A./(1+exp(-beta*F));
end

%%% stochastic choice
function [b] = choice(P)
    pp = rand();
    if pp<P
        b = 1;
    else
        b = 0;
    end
end

function [NLL] = nLL(THETA, dth, dcp, dc, Basis, lambda)
    
    %%% regularization
    if nargin < 6
        lambda = 0;
    end

    %%% Assume we parameterize in such way first
    alpha_h = THETA(1:4);
    alpha_dc = THETA(5:8);
    alpha_dcp = THETA(9:12);
    kappa_turn = THETA(13)^0.5;
    kappa_wv = THETA(14)^0.5;
    gamma = THETA(15);
    beta = 1;
    A = THETA(16);
    
    %%% kernel with basis
    K_h = (alpha_h*Basis');  % dth kernel
    K_dc = (alpha_dc*Basis');  % dC kernel
    K_dcp = (alpha_dcp*Basis');  % dCp kernel
    
    %%% turning decision
    d2r = 1;%pi/180;

    filt_dth = conv_kernel(abs(dth(1:end-1))*d2r,K_h);
    filt_dc = conv_kernel(dc(2:end),K_dc);
    P = NL(filt_dth + filt_dc, A);
    
    %%% weathervaning part
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    filt_dcp = conv_kernel(dcp(2:end),K_dcp);
    VM = C * exp(kappa_wv^2*cos(( filt_dcp - dth(2:end) )*d2r));  %von Mises distribution
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
    VM_turn = gamma*1/(2*pi) + (1-gamma)*VM_turn;  %%% revisit mixture inference !!!
    
    marginalP = (1-P).*VM + VM_turn.*P;
    
    NLL = -nansum(log(marginalP + 0*1e-10)) + lambda*sum(K_dc.^2);  % adding slope l2 regularization
end
