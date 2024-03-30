function [dth_pred, dr_pred] = model_predit_PAW(model_mmhat, xx, yy, gams_)
%%% given the fitted model, input xx,yy variables, and gams_ posterior states, 
%%% make predictions of the same length of time series for headings

lt = length(xx);
basis = model_mmhat.basis;
wind = length(basis);
xv = 1:wind;
dth_pred = zeros(1,lt);  % model prediction
dr_pred = zeros(1,lt);

for tt = length(xv)+1:lt
    %%% inputs
    C_t = xx(1, tt-wind+1:tt);
    dcp_t = xx(2, tt-wind+1:tt);
    ang_t = yy(1, tt-wind:tt-1);
    
    %%% choose state!
%     state_t = bb(tt);
    prob = gams_(:,tt);
    if rand()<prob(1)
        state_t = 1;
    else
        state_t = 2;
    end
    
    %%% unpack parameters given state
    if ndims(model_mmhat.wts)==3
        x = squeeze(model_mmhat.wts(:,:,state_t));
    else
        x = model_mmhat.wts;
    end
    K_ = x(1); B_ = x(2:5); Amp = x(6); tau = x(7); Amp_h = x(8); tau_h = x(9); K2_ = x(10);  gamma = x(11); A_=x(12); C_=x(13); b_dc=x(18); b_dcp=x(19);
    ks_ = x(14:15);  thetas_ = x(16:17); phis_ = [0,0];
    K_h_rec = Amp_h*exp(-xv/tau_h);
    K_dc_rec = B_*basis';
    K_dcp_rec = Amp*exp(-xv/tau);
    
    %%% modeling
    proj_dc = fliplr(K_dc_rec)*C_t';
    proj_dth = fliplr(K_h_rec)*abs(ang_t)';
    dc_dth = proj_dc + proj_dth + b_dc;
    wv = ((fliplr(K_dcp_rec)*dcp_t') + b_dcp) + (vmrand(0,K_))*180/pi;
    Pturns = (A_-C_) ./ (1 + exp( -(dc_dth)) + 0) + C_; %+sb
    
    if rand < Pturns*7/14
        beta = 1;
        dr = gamrnd(ks_(1), thetas_(1), 1, 1);
    else
        beta = 0;
        dr = gamrnd(ks_(2), thetas_(2), 1, 1);
    end

    if rand < gamma
        rt = beta*(vmrand(pi,1*K2_)*180/pi);
    else
        rt = beta*(vmrand(0,0)*180/pi);
    end
    dth = wv+rt;
    dth = wrapToPi(dth*pi/180)*180/pi;
            
    dth_pred(tt) = dth;
    dr_pred(tt) = dr;
end

dth_pred = dth_pred(wind:end);  % remove window
dr_pred = dr_pred(wind:end);

end