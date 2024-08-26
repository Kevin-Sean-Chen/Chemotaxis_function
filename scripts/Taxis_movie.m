%%Taxis_movie

%% load data (old)
% track = Tracks(14).Path;
% wim = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20220203_GWN_app+_sparse_11mM_9_30ml_400air/Data20220203_144516/individual_worm_imgs/worm_14.mat');
% wim = wim.worm_images;

%%
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low.mat');
Fcon = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/OdorFx_low.mat');
Fcon = Fcon.F;
M = Cmap.vq1;
M = fliplr(flipud(M));

%% load data (new)
% folder_name = '/projects/LEIFER/Kevin/Data_learn/N2/20231215_GWN_app_N2_11mM_35_400/Data20231215_132956';
folder_name = '/projects/LEIFER/Kevin/Data_learn/N2/20230721_GWN_app_N2_11mM_35_400/Data20230721_145443';
% for ii = 101:200
worm_id = 114;  %%33, 68
temp = load([folder_name,'/analysis/SmoothX.mat']);
x_pos = temp.values{worm_id};
temp = load([folder_name,'/analysis/SmoothY.mat']);
y_pos = temp.values{worm_id};
figure;
imagesc(M); hold on
plot(x_pos, y_pos, 'k')
% pause();
% end

wim = load([folder_name,'/individual_worm_imgs/worm_',num2str(worm_id),'.mat']);
wim = wim.worm_images;

%% full view version
time_window = [1:975];%[500:1500];
space_x = [300: 600]; %[950: 1250];
space_y = [1600:1900];  %[700: 1100];
range_xy = -35:34;


figure;

vidfile = VideoWriter('chemotaxis_movie2.avi');
open(vidfile);

original_cmap = colormap;  % Save the original colormap
colormap([1, 1, 1; original_cmap]);  % Add white for NaN
for tt = 1:3:length(time_window)
    
    ti = time_window(tt);
    
    %%% ego-centric worm
    worm_it = fliplr(im2double(wim(:,:,ti)));
    worm_it(find(worm_it>0)) = NaN;
    worm_it = fliplr(worm_it);
    
    %%% add to large scale
    x_it = round(x_pos(ti));
    y_it = round(y_pos(ti));
    M_temp = M;
    M_temp(y_it+range_xy, x_it+range_xy) = M_temp(y_it+range_xy, x_it+range_xy) + worm_it;
    M_temp = M_temp(space_y, space_x);
    
    %%% plot
    imagesc(M_temp(1:end-0,:));% hold on;
    set(gca, 'YDir', 'reverse');
    xticks([]); yticks([]); set(gcf,'color','w');
    colorbar('Fontsize',15);
    
%     colormap([1,1,1; colormap])
%     caxis([min(M_temp(:)),  max(M_temp(:))])
%     set(gca, 'CLim', [min(M_temp(:)) max(M_temp(:))]);
%     pause();
    
    drawnow
    F(tt) = getframe(gcf); 
    writeVideo(vidfile,F(tt));
    
end
close(vidfile)

%% ego-centric version
time_window = [500:1500];  %[900:1500];
range_xy = -35:34;

figure;
% vidfile = VideoWriter('chemotaxis_movie.avi');
% open(vidfile);

for tt = 1:7:length(time_window)
    
    %%% concentration
    ti = time_window(tt);
%     xyt = track(2,:);
%     center_c = Fcon(xyt(1), xyt(2));
%     range_x = range_xy + xyt(1);
%     range_y = range_xy + xyt(2);
%     [Xq,Yq] = meshgrid(range_x, range_x);
%     vi = Fcon(Xq,Yq);
%     vi = (flipud((vi')));

    %%% background sensory
    vi = M(round(x_pos(ti))+range_xy, round(y_pos(ti))+range_xy);
    %%% worm
    worm_it = fliplr(im2double(wim(:,:,ti)));
    worm_it(find(worm_it>0)) = NaN;
    
    
    %%% plot
    imagesc(vi+worm_it);% hold on;
    xticks([]); yticks([]); set(gcf,'color','w');
    colorbar('Fontsize',15);
%     caxis([21.5, 22.5]);
    pause();
%     drawnow
%     F(tt) = getframe(gcf); 
%     writeVideo(vidfile,F(tt));
    
end
% close(vidfile)
%%
bin = 1;
poly_degree = 3;
filt = 5;
x_smooth = smooth(track(:,1), filt,'sgolay',poly_degree);
y_smooth = smooth(track(:,2), filt,'sgolay',poly_degree);
temp = [x_smooth'; y_smooth']';
subs = temp(1:bin:end,:);
vecs = diff(subs);
Ct = zeros(1,length(vecs));
dtht = Ct*1;
for pp = 2:length(vecs)
    Ct(pp) = Fcon(subs(pp,1), subs(pp,2));
    dtht(pp) = angles(vecs(pp-1,:)/norm(vecs(pp-1,:)),vecs(pp,:)/norm(vecs(pp,:)));
end
time_ = [1:length(vecs)]/14*bin;
figure;
subplot(211); plot(Ct(2:end-1)); ylabel('ppm');
subplot(212); plot(dtht(2:end-1)); ylabel('d\theta'); xlabel('time (1/14 s)')

%% visualization
time_window = [900:1500];
range_xy = -35:34;

figure;
for tt = 1:length(time_window)
    
    %%% concentration
    ti = time_window(tt);
    xyt = track(2,:);
    center_c = Fcon(xyt(1), xyt(2));
    range_x = range_xy + xyt(1);
    range_y = range_xy + xyt(2);
    [Xq,Yq] = meshgrid(range_x, range_x);
    vi = Fcon(Xq,Yq);
    vi = (flipud((vi')));
    
    %%% worm
    worm_it = fliplr(im2double(wim(:,:,ti)));
    worm_it(find(worm_it>0)) = NaN;
    
    
    %%% plot
    imagesc(vi+worm_it);% hold on;
    xticks([]); yticks([]); set(gcf,'color','w');
    colorbar('Fontsize',15);
    caxis([21.5, 22.5]);
    pause();
    
end

%% 
figure
tt = 500
ti = time_window(tt);
xyt = track(2,:);
center_c = Fcon(xyt(1), xyt(2));
range_x = range_xy + xyt(1);
range_y = range_xy + xyt(2);
[Xq,Yq] = meshgrid(range_x, range_x);
vi = Fcon(Xq,Yq);
vi = (flipud((vi')));

%%% worm
worm_it = fliplr(im2double(wim(:,:,ti)));
worm_it(find(worm_it>0)) = NaN;


%%% plot
imm = flipud(vi)+fliplr(worm_it);
h=imagesc(imm);% hold on;
set(h, 'AlphaData', ~isnan(imm));set(gca,'color',0*[1 1 1]);
