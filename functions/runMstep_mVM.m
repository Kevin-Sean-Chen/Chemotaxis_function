function mm = runMstep_mVM(mm, xx, yy, gams, mask)
% mm = runMstep_LinGauss(mm,xx,yy,gams)
%
% Run m-step updates for Gaussian observation model
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%        .wts  [1 K] - per-state slopes
%        .vars [1 K] - per-state variances
%    xx [d T] - input (design matrix)
%    yy [1 T] - outputs
%  gams [K T] - log marginal sate probabilities given data
%
% Output
% ------
%  mmnew - new model struct

% normalize the gammas to sum to 1
gamnrm = gams./(sum(gams,2)+1e-8);  
nStates = size(gams,1);

%%% loading parameters
Basis = mm.basis;
lambda = mm.lambda;
dth = yy;
dc = xx(1,:);
dcp = xx(2,:);
nB = size(Basis,2);

% setup fmincon
% lfun = @(x)nll_mVM(x, dth, dcp, dc, gams, Basis, lambda, mask);
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','iter');
% opts.Algorithm = 'sqp';
LB = [1e-0, ones(1,nB)*-inf, -inf, 1e-0, -inf, 1e-1, 1e-0, 0.1];
UB = [30, ones(1,nB)*inf, inf, 100, inf, 100, 20, 1];
% prs0 = rand(1,10);
prs0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.] ;
% prs0 = [9.9763  -0.5343  -0.0776   0.1238  -0.0529   0.5335   7.7254  367.3817  0.1990  1.0000  0.1000]; %single mGLM
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);

for jj = 1:nStates
    lfun = @(x)nll_mVM(x, dth, dcp, dc, gamnrm(jj,:), Basis, lambda, mask);
%     prs0 = prs0 + prs0.*randn(1,length(UB))*0.5;
    prs0 = mm.wts(:,:,jj); %+ mm.wts(:,:,jj).*(2*(rand(1,length(UB))-0.5))*0.5;  %from last time!
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
%     x = fminunc(lfun, prs0, opts);
    mm.wts(:,:,jj) = x; % weighted optimization
end


end

function [nll] = nll_mVM(THETA, dth, dcp, dc, gams, Basis, lambda, mask)
    
    %%% Assume we parameterize in such way first
    kappa_wv = THETA(1)^0.5;      % variance of von Mises
    alpha_dc = THETA(2:5);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    Amp_dcp = THETA(6);       % amplitude for kernel for dcp normal concentration difference
    tau_dcp = THETA(7);       % time scale for dcp kernel
    Amp_h = THETA(8);         % amplitude of turning history kernel
    tau_h = THETA(9);         % time scale for dth history kernel
    kappa_turn = THETA(10)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(11);        % weight for uniform angle in the turn

    %%% turning decision
    %[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(5, [0, 10], 1.5);
    K_dc = alpha_dc * Basis';  %reconstruct with basis
    K_win = 1:length(K_dc);
    h_win = 1:length(K_dc)*1;
    K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel (make longer~)
    filt_dth = conv_kernel(abs(dth), K_h);
    filt_dc = conv_kernel(dc, K_dc);
    P = 1 ./ (1 + exp( -(filt_dc + filt_dth + 0))) +0;  %sigmoid(A_,B_,dc); 

    %%% weathervaning part
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    K_dcp = Amp_dcp * exp(-K_win/tau_dcp);    % dcp kernel
    filt_dcp = conv_kernel(dcp, K_dcp);
    VM = C * exp(kappa_wv^2*cos(( dth -filt_dcp )*d2r));  %von Mises distribution

    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth*d2r - pi)));  %test for non-uniform turns (sharp turns)
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;

    %%% marginal probability
    marginalP = (1-P).*VM + VM_turn.*P;
    nll = -nansum( mask .* gams .*( log( marginalP + 1*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));% + 0.1*sum((
end
