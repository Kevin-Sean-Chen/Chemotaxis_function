%% worm_FDT
%%%
% This script can be used after we have the 'Data' structure for chemotaxis
% analysis, which is a track-based structure with fields such as xy, dth,
% dc, dcp and others.
% The idea is to analyze these tracks with simple stat-mech methods and
% test for fluctuation-dissipation relations if possible
% If we need to exclude turning behavior, the fitted parameters from the
% model are needed.
%%%

%% MASD scaling
tau = 200;
ld = length(Data);
masd = zeros(1,tau);

for d = 1:1
    dthi = Data(d).theta; %dth;
    
    %%% conditioning on smooth runs
    ang_fit = Data(d).dth;
    ddc_fit = Data(d).dc;
    filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
    dc_dth = filt_ddc + filt_dth;
    Pturns = A_ ./ (1 + exp( -(dc_dth + base_dc))) + C_;
    run_pos = find(Pturns<0.99);
    mask = ones(1,length(dthi))*nan;
    mask(run_pos) = 1;
    dthi = dthi.*mask;
    
    %%% compute MSAD
    for ti = 1:length(dthi)-tau
        for tj = 1:tau
            if isnan(dthi(ti+tj))~=1 && isnan(dthi(ti))~=1  %%%abs(Data(d).dth(ti))<60  && abs(Data(d).dth(tj))<60
%           if abs(dthi(ti+tj) - dthi(ti))>60
%               masd(tj) = masd(tj)+ (dthi(ti+tj) - dthi(ti))^2;%+0;
%           else
            masd(tj) = masd(tj) + (dthi(ti+tj) - dthi(ti))^2;
          end
        end
    end
    masd = masd/length(dthi);
end
masd = masd/ld;

figure;
plot([1:tau]*5/14, masd, 'k', 'LineWidth',2)
xlabel('\tau (s)')
ylabel('<(\theta_t - \theta_0)^2>')
set(gcf,'color','w'); set(gca,'Fontsize',20);

%% FDT test
dcps = [];  % perturbation through c^perp
rts = [];  % recover time with angle
wind_dcp = 100;  % prior window to calculate dcp fluctuation
wind_rt = 100;  % causal window into the future to look for response time\
tol_dth = 5;  % the tolerance for dth around zero

for d = 1:100 %ld
    dthi = Data(d).theta; %dth;
    angi = Data(d).dth;
    dcpi = Data(d).dcp;
    
    %%% conditioning on smooth runs
    ang_fit = Data(d).dth;
    ddc_fit = Data(d).dc;
    filt_ddc = conv_kernel(ddc_fit, K_dc_rec);
    filt_dth = conv_kernel(abs(ang_fit), K_h_rec);
    dc_dth = filt_ddc + filt_dth;
    Pturns = A_ ./ (1 + exp( -(dc_dth + base_dc))) + C_;
    run_pos = find(Pturns<0.99);
    mask = ones(1,length(dthi))*nan;
    mask(run_pos) = 1;
    dthi = dthi.*mask;
    dcpi = dcpi.*mask;
    
    for ti = wind_dcp+1:1:length(dthi)-wind_rt
        dcps = [dcps std(dcpi(ti-wind_dcp:ti))];
%         [aa,bb] = min(find(abs(angi(ti:ti+wind_rt))<tol_dth)); %find the first time-step that returns around zero
        temp = find(abs(angi(ti:ti+wind_rt))<tol_dth);
        rts = [rts temp(1)];
    end
    d
end

figure;
plot(dcps, rts*5/14, 'o')
xlabel('fluctuation <(\delta C^\perp)^2>')
ylabel('response time (s)')
set(gcf,'color','w'); set(gca,'Fontsize',20);

