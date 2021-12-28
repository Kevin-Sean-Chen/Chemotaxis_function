%%% Track with Concentration
clear
clc

%% load tracks and landscape
Cmap = load('Landscape.mat');
Fcon = load('OdorFx.mat');
Fcon = Fcon.F;
Cmap = Cmap.vq1;
Cmap = fliplr(flipud(Cmap));  %%% for inverted camper!

% Paths = load('/projects/LEIFER/Kevin/20211013_biased_110mM/Data20211013_171205/analysis/Path.mat');
Paths = load('/projects/LEIFER/Kevin/20211013_biased_110mM/Path_1028.mat');
Paths = Paths.values;
% Paths = load('/projects/LEIFER/Kevin/20211013_biased_110mM/Data20211013_171205/tracking_deleted_tracks.mat');
% Paths = Paths.deleted_tracks;
Ntrack = length(Paths);

%% visualize
figure()
imagesc(Cmap)
hold on
for ii = 1:1:Ntrack
    plot(Paths{ii}(:,1), Paths{ii}(:,2), 'k')
%     plot(Paths(ii).Path(:,1), Paths(ii).Path(:,2), 'k')
end

%% without centerlines
del_cls = load('/projects/LEIFER/Kevin/20211029_GWN_app+_MEK110mM_40ml/Data20211029_175425/centerline_deleted_tracks.mat');
del_cls = del_cls.deleted_tracks;
figure()
imagesc(Cmap)
hold on
for ii = 1:length(del_cls)
    temp = del_cls(ii).Path;
    plot(temp(:,1), temp(:,2))
    hold on
end

%% make concentration profile for each track
Cons = {};
for ii = 1:Ntrack
    xy = round(Paths{ii});  %round for now... maybe extrapolate later
    Cs = Cmap(sub2ind(size(Cmap), xy(:,2), xy(:,1)));
%     Cs = Cmap([xy(:,1); xy(:,2)]');
    Cons{ii} = Cs;
end

%% plots
figure()
for ii = 1:Ntrack
    temp = Cons{ii};
%     plot(temp/max(temp))
    plot(temp - temp(1))
    hold on
end

%% %%% Test with angles %%%
Vels = load('/projects/LEIFER/Kevin/20211013_biased_110mM/Data20211013_171205/analysis/Velocity.mat');
Vels = Vels.values;
figure()
Ts = [];
Cs = [];
jj = 1;
for ii = 1:Ntrack
    temp = Cons{ii}-Cons{ii}(1);
    plot(temp, Vels{ii})
    hold on
    
    %%% record turning and concentration
    pos = find(Vels{ii}<0.01);  %reversal
%     if isempty(pos)~=1
    tur = zeros(1,length(temp));
    tur(pos) = 1;
    Ts = [ Ts, tur ];
    Cs = [ Cs, temp'];
end

%% fitting
dCs = diff(Cs);
[b,dev,stats] = glmfit(dCs, Ts(1:end-1)', 'binomial', 'link', 'probit')
xx = linspace(min(dCs), max(dCs), 5);
yfit = glmval(b,xx,'logit');
plot(xx,yfit,'-')
hold on
plot(dCs,Ts(1:end-1),'ro')
