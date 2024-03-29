function [dth] = P_dth_dC(prs,dcp_samp,ddc_samp)
    
    TT = length(dcp_samp);
    dth = zeros(1,TT);
    
    % unwrap parameters
    K_ = prs(1); 
    A_ = prs(2); 
    B_ = prs(3:6); 
    C_ = prs(7); 
    Amp = prs(8); 
    tau = prs(9); 
    K2_ = prs(10);
    
    %construction
    nB = 4;
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 10], 1.3);
    kappa = (1/K_)^0.5*(180/pi);  
    A = A_*1;
    xx = 1:length(cosBasis);
    Kdcp = Amp*exp(-xx/tau);
    Kddc = B_ * cosBasis';
    Pturn_base = C_;
    kappa2 = (1/K2_)^0.5*(180/pi);  
    
    % dynamics
    P_t = A ./ (1 + exp( conv( ddc_samp, Kddc, 'same' ))) + Pturn_base;  %biased random walk
    wv = -conv(dcp_samp, Kdcp, 'same') + kappa*randn(1,TT);  %weathervaning
    
    
    pos = find(rand(1,TT)< P_t);
    betas = zeros(1,TT);
    betas(pos) = 1;
    
    rt = betas.*( (randn(1,TT)*kappa2-180)*1. + 0.*(rand(1,TT)*360-180) );
    dth = wv + rt;
    
    dth = (mod(dth.*pi/180+pi,pi*2)-pi)./pi*180;
%     pos_l = find(dth>180);
%     pos_s = find(dth<180);
%     dth(pos_l) = dth(pos_l)-180;
%     dth(pos_s) = dth(pos_s)+360;
%     if dth>180; dth = dth-180; end;  if dth<-180; dth = dth+360; end  %within -180~180 degree range
    
end