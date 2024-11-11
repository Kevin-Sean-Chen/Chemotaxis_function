% bearing_samp
%%% sampling the bearing angle to compute alignment and errors in histrogram

%% load stored data
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare/bearing_data.mat');
temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/bearing_data.mat');
bt_data = temp.Bpairs(2,:);
bpair_data = temp.Bpairs;
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare/bearing_staPAW2.mat');
temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/bearing_staPAW.mat');
bt_stapaw = temp.Bpairs(2,:);
bpair_stapaw = temp.Bpairs;
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare/bearing_woF.mat');
temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/bearing_woF.mat');
bt_wof = temp.Bpairs(2,:);
bpair_wof = temp.Bpairs;
% temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare/bearing_dpaw.mat');
temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/bearing_dPAW.mat');
bt_dpaw = temp.Bpairs(2,:);
bpair_dpaw = temp.Bpairs;

temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/bearing_turn_in_state.mat');
bt_turnin = temp.Bpairs(2,:);
bpair_turnin = temp.Bpairs;
%%
temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/turns_staPAW2.mat');
tr_stapaw = temp.Bpairs(2,:);
tpair_stapaw = temp.Bpairs;
temp = load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/bearing/model_compare_0513/turns_data.mat');
tr_data = temp.Bpairs(2,:);
tpair_data = temp.Bpairs;

%% sampling for alignment
rng(42);
nbin = 12;  % number of bins
nsamp = 2400/2;  % number of samples
nrep = 50;  % repeats
b_samp = zeros(4, nrep);

for rr = 1:nrep
    temp = bt_data(randperm(length(bt_data)));
    b_samp(1,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    temp = bt_stapaw(randperm(length(bt_stapaw)));
    b_samp(2,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    temp = bt_wof(randperm(length(bt_wof)));
    b_samp(3,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    temp = bt_dpaw(randperm(length(bt_dpaw)));
    b_samp(4,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    
end

%% figure;
figure;
plot(repmat([1:4]',1,nrep)+randn(4,nrep)*0.09, b_samp,'k.','MarkerSize',13)
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('P(B aligned)'); xlabel('models')
xticks([1:4]); xticklabels({'data', 'staPAW', 'HMM', 'dPAW'}); xlim([0.5,4.5])

%% sampling for alignment.... for all
rng(42);
nbin = 12;  % number of bins
nsamp = 2400/2;  % number of samples
nrep = 50;  % repeats
b_samp = zeros(6, nrep);

for rr = 1:nrep
    temp = bt_data(randperm(length(bt_data)));
    b_samp(1,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    temp = tr_data(randperm(length(tr_data)));
    b_samp(2,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    
    temp = bt_stapaw(randperm(length(bt_stapaw)));
    b_samp(3,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    temp = tr_stapaw(randperm(length(tr_stapaw)));
    b_samp(4,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    
    temp = bt_wof(randperm(length(bt_wof)));
    b_samp(5,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    temp = bt_dpaw(randperm(length(bt_dpaw)));
    b_samp(6,rr) = length(find(abs(temp(1:nsamp))<90))/nsamp;
    
end

%% figure;
figure;
plot(repmat([1:6]',1,nrep)+randn(6,nrep)*0.09, b_samp,'k.','MarkerSize',13)
set(gcf,'color','w'); set(gca,'Fontsize',20); ylabel('P(B aligned)'); xlabel('models')
xticks([1:6]); xticklabels({'data', 'turn', 'staPAW','turn', 'HMM', 'dPAW'}); xlim([0.5,6.5])

%% histograms with uncertainty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test one model/data at a time
this_b = bt_data; this_pair = bpair_data;  % data
this_b = bt_stapaw; this_pair = bpair_stapaw;  % staPAW
this_b = bt_wof; this_pair = bpair_wof;  % woF
this_b = bt_dpaw; this_pair = bpair_dpaw;  % dPAW

this_b = tr_data; this_pair = tpair_data;  % turn for data
this_b = tr_stapaw; this_pair = tpair_stapaw;  % turn for staPAW

this_b = bt_turnin; this_pair = bpair_turnin;
% this_pair = Bpairs; this_b = Bpairs(2,:);  % external data, for mutants and others

rng(42)

%%
nbin = 12;  % number of bins
nsamp = 2400/2;%1680/1;  % number of samples (10min)
nrep = 100;  % repeats
b_hist = zeros(nrep, nbin);
samp_record = [];
for rr = 1:nrep
    temp = this_b(randperm(length(this_b)));
    samp_b = temp(1:nsamp);
    [pb,midb] = bear2hist(samp_b, nbin);
    b_hist(rr,:) = pb;
    samp_record = [samp_record samp_b];
end

test = samp_record;
length(find(abs(test)<90))/length(test) 
%%
figure
% bar(midb, mean(b_hist),'k','facealpha', 0.7); hold on
% errorbar(midb, mean(b_hist), std(b_hist)/sqrt(1));
histogram(samp_record, nbin, 'FaceColor', 'b', 'EdgeColor', 'none', 'Normalization', 'probability'); hold on;
fill([midb, fliplr(midb)], [mean(b_hist)-std(b_hist), fliplr(mean(b_hist)+std(b_hist))], 'k','edgecolor', 'none', 'facealpha', 0.4);  
%%% control
bc = bear2control(this_pair, 100);
[pb,midb] = bear2hist(bc, nbin);
plot(midb, pb/sum(pb),'k', 'LineWidth',3)
%%% work on filled area + histbar + shuffled line...
xlim([-180, 180]); ylim([0, .15]);
set(gcf,'color','w'); set(gca,'Fontsize',20); xticks([])

%% function
function [pb,midb] = bear2hist(B, nbin)
    H = histogram(B, nbin, 'Normalization', 'pdf','Visible', 'off');
    binEdges = H.BinEdges;  
    pb = H.BinCounts;
    pb = pb/sum(pb);  
    midb = (binEdges(2:end)+binEdges(1:end-1))/2;
%     close
end

function [Bc] = bear2control(Bpairs, shuffle_rep)
    Bc  = [];
    dB = -diff(Bpairs,1);
    for ss = 1:shuffle_rep
        Bc = [Bc dB(randperm(length(dB))) + Bpairs(1,:)];
    end
    Bc = wrapToPi(Bc*pi/180)*180/pi;
end