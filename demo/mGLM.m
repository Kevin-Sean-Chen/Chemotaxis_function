% mGLM
%%% mixture GLM debugging
%% define variables
lt = 10000; % length of simulation / data
nB = 4;  % number of basis function for the kernel
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3); % basis function
%%% true params
beta = 2; % nonlinear parameter
alpha_h = [-4:-1]*.1;  % angle history kernel coefficient
alpha_dc = [1,4,-2,-1]*1;  % dC kernel coefficient
alpha_dcp = [-4:-1]*.01;  % dCp kernel coefficient
base = 0;  %baseline
kappa_turn = 5;  % turning angle variance
kappa_wv = 10;  % weather-vaning angle variance

%% generate data
K_h = fliplr(alpha_h*cosBasis');  % dth kernel
K_dc = fliplr(alpha_dc*cosBasis');  % dC kernel
K_dcp = fliplr(alpha_dcp*cosBasis');  % dCp kernel
lc = 5;  %length of smoothing
dC = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dC stimulus vector
dCp = conv(randn(1,lt),ones(1,lc),'same')/lc;  % dCp stimulus vector
dth = zeros(1,lt);
turns = zeros(1,lt);
F = dth*0;
pad = length(K_s);
for tt=pad:lt
    F(tt) = dC(tt-pad+1:tt)*K_dc' + abs(dth(tt-pad+1:tt))*K_h';  % linear filtering
    turns(tt) = choice(NL(F(tt)+base,beta));  % nonlinearity and binary choice
    dth(tt) = turns(tt)*circ_vmrnd(pi,kappa_turn,1) + (1-turns(tt))*circ_vmrnd(dCp(tt-pad+1:tt)*K_dcp',kappa_wv,1);  % angle drawn from mixture of von Mesis
end

%% MLE inference
lfun = @(x)nLL(x, dth, dCp, dC, cosBasis);  % objective function
opts = optimset('display','iter');
num_par = 14;
LB = [ones(1,12)*-10, 0, 0]*1;
UB = [ones(1,12)*10, 20, 20]*1;
% prs0 = [alpha_h, alpha_dc, alpha_dcp, kappa_turn, kappa_wv];%
prs0 = rand(1,num_par);
[x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);  % constrained optimization
% [x,FVAL,EXITFLAG,OUTPUT,GRAD,HESSIAN] = fminunc(lfun, prs0, opts);

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

%% functions
%%% nonlinear function
function [P] = NL(F,beta)
    P = 1./(1+exp(-beta*F));
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
    beta = 2;
    
    %%% kernel with basis
    K_h = (alpha_h*Basis');  % dth kernel
    K_dc = (alpha_dc*Basis');  % dC kernel
    K_dcp = (alpha_dcp*Basis');  % dCp kernel
    
    %%% turning decision
    d2r = 1;%pi/180;
%     padding = ones(1,floor(length(K_h)/2));
%     filt_dth = conv([padding, abs(dth)*d2r], K_h, 'same');
%     filt_dth = filt_dth(1:length(dth));
    filt_dth = conv_kernel(abs(dth)*d2r,K_h);
    filt_dc = conv_kernel(dc,K_dc);
    P = NL(filt_dth + filt_dc, beta);
%     P = 1./(1 + exp(-beta*(filt_dth + filt_dc)));
    
    %%% weathervaning part
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    filt_dcp = conv_kernel(dcp,K_dcp);
    VM = C * exp(kappa_wv^2*cos(( filt_dcp - dth )*d2r));  %von Mises distribution
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth*d2r - pi)));  %test for non-uniform turns (sharp turns)
%     gamma = 1;
%     VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    
    marginalP = (1-P).*VM + VM_turn.*P;
    
%     lambda = 10;
    NLL = -nansum(log(marginalP + 0*1e-10));  % adding slope l2 regularization
end
