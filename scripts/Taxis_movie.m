%%Taxis_movie

%% load data
track = Tracks(14).Path;
wim = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/20220203_GWN_app+_sparse_11mM_9_30ml_400air/Data20220203_144516/individual_worm_imgs/worm_14.mat');
wim = wim.worm_images;
Cmap = load('/home/kschen/github/OdorSensorArray/OSA_MFC_PID_scripts/Landscape_low.mat');
Fcon = load('/home/kschen/github/OdorSensorArray/OSA_MFC_PID_scripts/OdorFx_low.mat');
Fcon = Fcon.F;
M = Cmap.vq1;
M = (flipud(M));

%%
time_window = [900:1500];
range_xy = -35:34;

figure;
vidfile = VideoWriter('chemotaxis_movie.avi');
open(vidfile);

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
%     pause();
    drawnow
    F(tt) = getframe(gcf); 
    writeVideo(vidfile,F(tt));
    
end
close(vidfile)
%%
bin = 1;
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