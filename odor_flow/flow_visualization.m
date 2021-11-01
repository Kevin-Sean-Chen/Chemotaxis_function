%%%flow_visualization

%% through time
% dir_ = '/tigress/LEIFER/Kevin/20190416/Data20190416_174432/';
% dir_ = '/tigress/LEIFER/Kevin/20190416/Data20190416_181432/';
dir_ = '/tigress/LEIFER/Kevin/20190416/Data20190418_104052/';
% dir_ = '/tigress/LEIFER/Kevin/20190419/Data20190419_134331/';

%%%with flow meter
% dir_ = '/tigress/LEIFER/Kevin/20190502/Data20190502_151721/';
% dir_ = '/tigress/LEIFER/Kevin/20190507/Data20190507_112612/';

%%%with flow meter and vac
% dir_ = '/tigress/LEIFER/Kevin/20190706/Data20190706_132342/';
list = dir(dir_);
for ii=3:50:2500
    M = imread([dir_,list(ii).name]); %raw unit8 images
    imagesc(M); 
    title(ii);
    drawnow; 
    pause();
end

%% difference
allback = im2double(imread([dir_,list(5).name])); 
for b=6:11
    temp = im2double(imread([dir_,list(b).name]));
    allback = cat(3,allback,temp); 
end
back = sum(allback,3)/size(allback,3);
% back = im2double(imread([dir_,list(200).name]));  %initial frame for backgrounbd subtraction
for ii=3:20:2500
    M = im2double(imread([dir_,list(ii).name])); 
    imagesc(M-back); 
    title(ii);
    drawnow; 
    pause();
end

%% stability

%% spatial profile
ii = 1100;  %1700
M = im2double(imread([dir_,list(ii).name]));
S = M-back;
c = 1;
col = hsv(length(1000:100:2000));
for cc = 1000:100:2000
    sc = S(:,cc);
    sc2 = im2double(sc);
    plot(smooth(sc2,100,'moving'),'color',col(c,:))
    hold on
    c = c+1;
end

%% scan spatial profile
pix2mm = 1;%/31;  %pixel number to mm measured from the ruler image

for ii = 5:10:length(list)

M = im2double(imread([dir_,list(ii).name]));
subplot(1,3,1); imagesc(M-back,'XData',[0 size(M,2)*pix2mm],'YData',[0 size(M,1)*pix2mm]); title(ii);
S = M-back;
S(S<0) = 0;  %rectify
c = 1;
col = hsv(length(1000:100:2000));
for cc = 1000:100:2000
    sc = S(:,cc);
    sc2 = im2double(sc);
    subplot(1,3,2); plot([1:length(sc2)]*pix2mm,smooth(sc2,100,'moving'),'color',col(c,:))
    hold on
    c = c+1;
end
hold off
sc = S(900,:);  %900 %center-line
sc2 = im2double(sc);
subplot(1,3,3); plot([1:length(sc2)]*pix2mm,smooth(sc2,100,'moving'))
pause();

end
