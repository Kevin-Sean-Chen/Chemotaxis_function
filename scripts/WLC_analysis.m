%%% WLC_analysis
%%
% Paths = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/App+_tracks.mat');
Paths = load('/projects/LEIFER/Kevin/Data_learn/N2/data_analysis/Tracks_nai.mat');
Tracks = Paths.Tracks;
% Tracks = deleted_tracks;

%% load force (odor) field and compute gradients
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
M = Cmap.vq1;
M = fliplr(flipud(M)); 
n_smooth = 100;
[fx,fy] = gradient((conv2(M,ones(n_smooth,n_smooth),'same')/n_smooth^2),1);

%%
figure;
%%% pre-processing parameters
poly_degree = 3;  %polynomial fit for the tracks
bin = 7;  %temporal binning  %~0.5 s
filt = 7;  %filtering tracks
l_window = 2;  %lag time

%%% gather time series
dTs = {};  %angle change
dFs = {};  %effective forcing
dDs = {};  %displacement
dRs = {};  %diract calculation of extension

%%% "goal"
[maxValue, linearIndex] = max(M(:));
[max_y, max_x] = ind2sub(size(M), linearIndex);

for c = 1:length(Tracks) %length(Paths)
    
    %%% process path to angles and distances
%     temp = Paths{c};  %for loading deleted tracks
%     temp = Tracks(c).Path;  %for loading saved tracks

    runs = Tracks(c).Runs;
    paths = Tracks(c).Path;
    
    %%% loop through runs
    for rr = 1%1:size(runs,1)
%         temp1 = paths(runs(rr,1):runs(rr,2),:); 
      temp1 = paths*1;
      
%     temp1 = zeros(round(size(temp,1)/1),2);
%     temp1(:,1) = smooth(temp(1:round(length(temp)/1),1), filt,'sgolay',poly_degree);
%     temp1(:,2) = smooth(temp(1:round(length(temp)/1),2), filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    vecs = diff(subs);
    angs = zeros(1,size(vecs,1));    

    %%% for time step calculations
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    
    %%% global grad
    grad_i = [subs(end,1)-subs(1,1), subs(end,2)-subs(1,2)];
%     grad_i = [max_x-subs(1,1), max_y-subs(1,2)];
    
    %%% iterate through worms
    for dd = (l_window+1):length(angs)
%         grad_i = [fx(floor(subs(dd,2)),floor(subs(dd,1))), fy(floor(subs(dd,2)),floor(subs(dd,1)))];
        
        %%% angle function
%         angs(dd) = angles(vecs(dd-1,:)/norm(vecs(dd-1,:)),vecs(dd,:)/norm(vecs(dd,:)));
%         angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        %%% angle is between gradients!!
        %%%
        angs(dd) = angles(vecs(dd,:)/norm(vecs(dd,:)), grad_i/norm(grad_i));
        %%%
        
        %%% perpendicular concentration change
        perp_dir = [-vecs(dd-l_window,2), vecs(dd-l_window,1)];
        perp_dir = perp_dir/norm(perp_dir);
        dCp(dd) = Fcon(subs(dd-l_window,1)+perp_dir(1)*1., subs(dd-l_window,2)+perp_dir(2)*1)...
             - Fcon(subs(dd-l_window,1)-perp_dir(1)*1, subs(dd-l_window,2)-perp_dir(2)*1);
        
        %forward concentration change
        dCs(dd) = Fcon(subs(dd,1),subs(dd,2)) - Fcon(subs(dd-l_window,1),subs(dd-l_window,2));
        
        %%%check displacement
        dds(dd) = norm(subs(dd,1)-subs(dd-l_window,1), subs(dd,2)-subs(dd-l_window,2));
        %norm(vecs(dd,:));

    end
    
    DC = (Fcon(temp1(end,2), temp1(end,1)) - Fcon(temp1(1,2), temp1(1,1)));% / Fcon(temp1(1,2), temp1(1,1));
    Ri = sqrt((temp1(end,1)-temp1(1,1))^2 + (temp1(end,2)-temp1(1,2))^2);
    
    %remove zeros
    dCp = dCp((l_window+1):end);
    dCs = dCs((l_window+1):end);
    dds = dds((l_window+1):end);
    angs = angs((l_window+1):end);
    
    end
    
    dTs{c} = angs;
%     dFs{c} = dCs+dCp*0;  %orthogonal force decomposition!
    dFs{c} = DC;
    dDs{c} = dds;
    dRs{c} = Ri;
    
end

%% F-R plot
figure();
fss = zeros(1,length(Tracks));
rss = zeros(1,length(Tracks));
for ii = 1:length(Tracks)
    
    R = 0;
    L = 0;
    f = 0;
    angs = dTs{ii};
    dis = dDs{ii};
    dC = dFs{ii};
    for rr = 1:length(angs)
%         R = R + cos(angs(rr)*pi/180);  %extension R
%         L = L + dis(rr);  %length L
%         f = f + dC(rr);  %force f

        R = R + dis(rr)*cos(angs(rr)*pi/180);  %extension R

    end
    
    %%% new calculation
    f = dFs{ii};
    L = sum(dis);
%     R = dRs{ii};
    
    fss(ii) = f;
    rss(ii) = R;
    
    semilogx(abs(f),R/L,'ko')
%     plot(f,R/L,'o')
    hold on
end
xlabel('effective force \Sigma \DeltaC')
ylabel('extension <R>/L')

set(gcf,'color','w');
set(gca,'FontSize', 20);

%% 
figure()
f_ = 1:5000;%5*logspace(-5,3);
%1:5000;%abs(fss);  %model input
L = 800;%max(rss);  %maximum length
lp = 0.1;%persistant length for guess...
ll = L/lp;

%%% WLC prediction via stat-mech
% R_ = L*( 1 - .1/4./sqrt(f_*lp) .* coth(ll.*sqrt(f_*lp)) + 1/4./(f_*lp) .* 1/ll );
R_ = L*( 1-1/4./sqrt(f_*lp) );

semilogx(f_,R_)
legend('WCL')
hold on
semilogx(abs(fss),rss,'o')

xlabel('effective force \Sigma \DeltaC')
ylabel('extension <R>')

