% opto_staPAW_analysis
%%%
% analyze the fitted staPAW in response to opto impulses
%%%
%% load
neuron = 'AWC';  %'RIM' 'AWC'
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/odor_opto_AWC_vars.mat')
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/odor_opto_RIM_vars.mat')
filename = ['/projects/LEIFER/Kevin/Publications/Chen_states_2024/odor_opto_',neuron,'_vars.mat'];
load(filename)

%% compile time seris
% with Data structure
[xxf, yyf, alltrials, time] = data2xy(Data);
alldis = extractfield(Data, 'dis');  % Data
allopto = extractfield(Data, 'opto');  % Data
xxf = [xxf; allopto];
yyf = [yyf; alldis];

% testing
wind_test = [1:500000]; %[100000:148982];%500000:length(allas)];%max(wind):length(allas);
offset = min(wind_test)-1;
yy = yyf(:,wind_test);
xx = xxf(:,wind_test);
mask = maskf(wind_test);
[logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xx,yy,mask);
[~, ~, ~, alltime] = data2xy(Data);

%%
pz_trig = 1;
pre_t = -5;
post_t = 15;
time_v = [pre_t:5/14:0   5/14:5/14:post_t];
acu_l = length(pre_t:5/14:0);  cua_l = length(5/14:5/14:post_t);
wl = length(time_v);
opto_triggs_z = [];
opto_triggs_dth = [];
opto_triggs_dr = [];
opto_t = xx(3,:);  %opto time series
trigg_t = find(diff(opto_t(1:end-wl))>1);  % onset
for trs = 1:length(trigg_t)
    pos = trigg_t(trs);
%     pos_r = randi([acu_l   length(opto_t)-cua_l]);  % random control!
    opto_triggs_z = [opto_triggs_z; gams_(pz_trig,pos-acu_l:pos+cua_l-1)];% - gams_(pz_trig,pos_r-acu_l:pos_r+cua_l-1)];
    opto_triggs_dth = [opto_triggs_dth; abs(yy(1,pos-acu_l:pos+cua_l-1))];% - abs(yy(1,pos_r-acu_l:pos_r+cua_l-1))];
    opto_triggs_dr = [opto_triggs_dr; yy(2,pos-acu_l:pos+cua_l-1)*1/30];% - yy(2,pos_r-acu_l:pos_r+cua_l-1)*1/30];
end

%%
sim_v = time_v*0;
sim_v(16:30) = 1;
figure;
subplot(4,1,1); plot(time_v, sim_v)
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('opto'); title([neuron,'::ChR2'])
subplot(4,1,2); plot(time_v, mean(opto_triggs_z)); hold on; 
plot(time_v, mean(opto_triggs_z)+std(opto_triggs_z)/sqrt(trs),'--'); plot(time_v, mean(opto_triggs_z)-std(opto_triggs_z)/sqrt(trs),'--')
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('P(z)')
subplot(4,1,3); plot(time_v, mean(opto_triggs_dth)); hold on; 
plot(time_v, mean(opto_triggs_dth)+std(opto_triggs_dth)/sqrt(trs),'--'); plot(time_v, mean(opto_triggs_dth)-std(opto_triggs_dth)/sqrt(trs),'--')
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('|d\theta|')
subplot(4,1,4); plot(time_v, mean(opto_triggs_dr)); hold on; 
plot(time_v, mean(opto_triggs_dr)+std(opto_triggs_dr)/sqrt(trs),'--'); plot(time_v, mean(opto_triggs_dr)-std(opto_triggs_dr)/sqrt(trs),'--')
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('dr'); xlabel('time (s)')

%% conditional trigger!
pre_z = 2;
opto_triggs_z_cond = [];
opto_triggs_dth_cond = [];
opto_triggs_dr_cond = [];
opto_t = xx(3,:);  %opto time series
trigg_t = find(diff(opto_t(1:end-wl))>1);  % onset >1 off <-1
for trs = 1:length(trigg_t)
    pos = trigg_t(trs);
%     pos = randi([acu_l   length(opto_t)-cua_l]);  % random control!
    [~,pre_stim_z] = max(sum(gams_(:,pos-acu_l:pos),2));
    if pre_stim_z==pre_z
        opto_triggs_z_cond = [opto_triggs_z_cond; gams_(pre_z, pos-acu_l:pos+cua_l-1)];% - gams_(pre_z, pos_r-acu_l:pos_r+cua_l-1)];
        opto_triggs_dth_cond = [opto_triggs_dth_cond; abs(yy(1,pos-acu_l:pos+cua_l-1))];% 
        opto_triggs_dr_cond = [opto_triggs_dr_cond; yy(2,pos-acu_l:pos+cua_l-1)*1/30];% 
    end
end

%%
sim_v = time_v*0;
sim_v(16:30) = 1;
figure;
subplot(4,1,1); plot(time_v, sim_v)
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('opto'); title([neuron,'::ChR2 conditioned on state',num2str(pre_z)])
subplot(4,1,2); plot(time_v, mean(opto_triggs_z_cond)); hold on; 
plot(time_v, mean(opto_triggs_z_cond)+std(opto_triggs_z_cond)/sqrt(trs),'--'); plot(time_v, mean(opto_triggs_z_cond)-std(opto_triggs_z_cond)/sqrt(trs),'--')
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('P(z)')
subplot(4,1,3); plot(time_v, mean(opto_triggs_dth_cond)); hold on; 
plot(time_v, mean(opto_triggs_dth_cond)+std(opto_triggs_dth_cond)/sqrt(trs),'--'); plot(time_v, mean(opto_triggs_dth_cond)-std(opto_triggs_dth_cond)/sqrt(trs),'--')
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('|d\theta|')
subplot(4,1,4); plot(time_v, mean(opto_triggs_dr_cond)); hold on; 
plot(time_v, mean(opto_triggs_dr_cond)+std(opto_triggs_dr_cond)/sqrt(trs),'--'); plot(time_v, mean(opto_triggs_dr_cond)-std(opto_triggs_dr_cond)/sqrt(trs),'--')
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('dr'); xlabel('time (s)')
