function mm = runMstep_staPAW(mm, xx, yy, gams, mask)
% mm = runMstep_LinGauss(mm,xx,yy,gams)
%
% Run m-step updates for staPAW observation model
%
% Inputs
% ------
%   mm [struct] - model structure with params 
%        .wts  [D K] - per-state parameters for dPAW
%        .A - transition baseline
%        .wts_state [K K d] - transition kernels
%        .basis, .lambda, .loglifun, .loglitrans...
%    xx [2 T] - input (concentration and perpendicular concentration)
%    yy [2 T] - outputs (angle and speed)
%  gams [K T] - log marginal sate probabilities given data
%  mask [1 T] - logic mask for likelihood function
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
dth = yy(1,:);
dv = yy(2,:);
dc = xx(1,:);
dcp = xx(2,:);
nB = size(Basis,2);

% setup fmincon
% lfun = @(x)nll_mVM(x, dth, dcp, dc, gams, Basis, lambda, mask);
% [x,fval] = fminunc(lfun,randn(1,10));  %random initiation
% [x,fval,exitflag,output,grad,hessian] = fminunc(lfun,[500, 0.0, randn(1,6), -1, 100]+randn(1,10)*0.);  %a closer to a reasonable value

opts = optimset('display','final');%'iter');
% opts.Algorithm = 'sqp';
LB = [1e-0, ones(1,nB)*-inf, -inf, 1e-0, -inf, 1e-1, 1e-0,   0.  0   0    0   0   0 0 -10 -10];
UB = [100,  ones(1,nB)*inf,   inf,  100,  inf, 100,   100,   1  1  0.5  100 100  50 50 10  10];
% prs0 = rand(1,10);
prs0 = [10, randn(1,nB)*10, 10, 25, 10, 25, 5, 1.   0.5 0  1 1 1 1 .1 .1 0 0] ;
% [x,fval] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);
nparam = length(LB);
ineq = zeros(nparam,nparam);
ineq(:,12) = -1;  ineq(:,13)=1;
ineq_r = zeros(1,nparam);
    
for jj = 1:nStates
    lfun = @(x)nll_staPAW(x, dth, dcp, dc, dv, gamnrm(jj,:), Basis, lambda(jj), mask);
%     prs0 = prs0 + prs0.*randn(1,length(UB))*0.5;
    prs0 = mm.wts(:,:,jj); %+ mm.wts(:,:,jj).*(2*(rand(1,length(UB))-0.5))*0.5;  %from last time!
    [x,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,ineq,ineq_r,[],[],LB,UB,[],opts);
%     x = fminunc(lfun, prs0, opts);
    mm.wts(:,:,jj) = x; % weighted optimization
end

end

function [nll] = nll_staPAW(THETA, dth, dcp, dc, dv, gams, Basis, lambda, mask)
    
    %%% Assume we parameterize in such way first
    kappa_wv = THETA(1)^0.5;      % variance of von Mises
    alpha_dc = THETA(2:5);    % kernel for dC transitional concentration difference (weights on kerenl basis)
    Amp_dcp = THETA(6);       % amplitude for kernel for dcp normal concentration difference
    tau_dcp = THETA(7);       % time scale for dcp kernel
    Amp_h = THETA(8);         % amplitude of turning history kernel
    tau_h = THETA(9);         % time scale for dth history kernel
    kappa_turn = THETA(10)^0.5;   % vairance of the sharp turn von Mises
    gamma = THETA(11);        % weight for uniform angle in the turn
    A_ = THETA(12);           % max turn probability
    B_ = THETA(13);           % basline turn probability
    gamm_shapes = THETA(14:15);    % shape parameters for gamma
    gamm_scales = THETA(16:17);    % scale parameters for gamma
%     gamm_AR = THETA(18:19);        % AR weight for gamma
%     base_dc = THETA(20);      % bias of turning
%     base_dcp = THETA(21);     % bias of curving
    gamm_AR = [0 0];
    base_dc = THETA(18);      % bias of turning
    base_dcp = THETA(19);     % bias of curving

    %%% turning decision
    %[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(5, [0, 10], 1.5);
    K_dc = alpha_dc * Basis';  %reconstruct with basis
    K_win = 1:length(K_dc);
    h_win = 1:length(K_dc)*1;
    K_h = Amp_h * exp(-h_win/tau_h);  % dth kernel (make longer~)
    filt_dth = conv_kernel(abs(dth(1:end-1)), K_h);
    filt_dc = conv_kernel(dc(2:end), K_dc);
    P = (A_-B_) ./ (1 + exp( -(filt_dc + filt_dth + base_dc))) + B_;  %sigmoid(A_,B_,dc); 
%     P = P*14/5;
    
    %%% weathervaning part
    d2r = pi/180;
    C = 1/(2*pi*besseli(0,kappa_wv^2));  % normalize for von Mises
    K_dcp = Amp_dcp * exp(-K_win/tau_dcp);    % dcp kernel
    filt_dcp = conv_kernel(dcp(2:end), K_dcp);
    VM = C * exp(kappa_wv^2*cos(wrapTo180( dth(2:end) - filt_dcp - base_dcp)*d2r));  %von Mises distribution
%     VM = C * exp(kappa_wv^2*cos(( dth(2:end) - filt_dcp - base_dcp)*d2r));
    gamm_wv = gampdf(dv(2:end), gamm_shapes(2) + gamm_AR(2)*dv(1:end-1), gamm_scales(2));
    
    %%% turning analge model
    VM_turn = 1/(2*pi*besseli(0,kappa_turn^2)) * exp(kappa_turn^2*cos((dth(2:end)*d2r - pi)));  %test for non-uniform turns (sharp turns)
    VM_turn = (1-gamma)*1/(2*pi) + gamma*VM_turn;
    gamm_pr = gampdf(dv(2:end), gamm_shapes(1) + gamm_AR(1)*dv(1:end-1), gamm_scales(1));
    
    %%% marginal probability
    marginalP = (1-P).*VM.*gamm_wv + P.*VM_turn.*gamm_pr;
    nll_dth = -nansum( mask(2:end) .* gams(2:end) .* ( log(marginalP + 0*1e-10) ) ) + lambda*(1*sum((K_dc - 0).^2));
%     pos = (~isinf(log(dv)));
%     ll_gam = (theta_gamm(1)-1)*nansum(gams(2:end).*log(dv(2:end)).*pos(2:end)) - nansum(gams(2:end).*dv(2:end))/theta_gamm(2) - nansum(gams(2:end).*theta_gamm(1)*log(theta_gamm(2)) + gams(2:end).*gammaln(theta_gamm(1)));
    nll = nll_dth;%- ll_gam;
end
