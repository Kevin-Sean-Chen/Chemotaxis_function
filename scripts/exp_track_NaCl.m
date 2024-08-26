% exp_track_NaCl
%%% demo pretty tracks for salt chemotaxis

%% load data
load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt100_50.mat')
% load('/projects/LEIFER/Kevin/Data_salt/data_analysis/Data_salt0_50_0513.mat')

%% for salt
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
M = Cmap.vq1;
M = fliplr(flipud(M));
[rows, cols] = size(M);
[x_, y_] = meshgrid(linspace(0, 50, cols), 1:rows);  % for   0 to 50 mM
% [x_, y_] = meshgrid(linspace(100, 50, cols), 1:rows);  % for 100 to 50 mM
gradient_x = x_ * 1;
M = (y_*0+1) .* gradient_x;  
%figure; imagesc(M)
pix2mm = 1/30;

%% examples
data_i = Data([57:77]); %50-100   %62-77
% data_i = Data(10:30); %0-50  %40:70; %50:80

%%% plotting
figure()
ax1 = axes;
imagesc(ax1,M,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]);
set(gcf,'color','w'); set(gca,'Fontsize',20); xticks([]); yticks([]);
colormap()
hold on
ax2 = axes;
for ii = 1:length(data_i)
    xx = data_i(ii).xy(1,:);
    yy = data_i(ii).xy(2,:);
    ll = length(data_i(ii).xy);
    gg = linspace(0,1,ll);
%     plot(xx, yy, 'Color', [grayLevel grayLevel grayLevel]);
    patch(ax2, [xx nan]*pix2mm,[yy nan]*pix2mm,[gg nan],[gg nan], 'edgecolor', 'interp','LineWidth',2); 
    hold on
    plot(ax2,xx(1)*pix2mm, yy(1)*pix2mm,'g.', 'MarkerSize',25)
    plot(ax2,xx(end)*pix2mm, yy(end)*pix2mm,'r.', 'MarkerSize',25)
end
set(gca, 'YDir','reverse')
set(gca, 'XDir', 'reverse')  %%% for 100-50
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
c = gray;
colormap(ax2,c)
colormap(ax1)
xlim([1,3000].*pix2mm); ylim([1,2500].*pix2mm)
