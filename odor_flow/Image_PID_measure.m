%%Image_PID_measure

%%%with flow meter and vac
dir_ = '/tigress/LEIFER/Kevin/20190607/Data20190607_153433/';
dir_ = '/tigress/LEIFER/Kevin/20190614/Data20190614_125818/';
dir_ = '/tigress/LEIFER/Kevin/20190618/Data20190618_113509/';  %calibration (missing PIDread...)
dir_ = '/tigress/LEIFER/Kevin/20190622/Data20190622_154522/';  %calibration
dir_ = '/tigress/LEIFER/Kevin/20190622/Data20190622_164113/';  %second calibration after 30min
dir_ = '/tigress/LEIFER/Kevin/20190622/Data20190622_155953/';  %stationarity check
dir_ = '/tigress/LEIFER/Kevin/20190623/Data20190623_165139/';
dir_ = '/tigress/LEIFER/Kevin/20190625/Data20190625_183624/';  %stationarity throught time
dir_ = '/tigress/LEIFER/Kevin/20190706/Data20190706_135201/';
dir_ = '/tigress/LEIFER/Kevin/20190718/Data20190718_112626/';  %SPGC3 calibration
dir_ = '/tigress/LEIFER/Kevin/20190719/Data20190719_142349/';  %EtOH calibration
dir_ = '/tigress/LEIFER/Kevin/20190722/Data20190722_132808/';  %compare SGPC3 and SGP30 in flow chamber
dir_ = '/tigress/LEIFER/Kevin/20190812/Data20190812_152806/';  %for stability check with divider-like setup and in the new chamber
list = dir(dir_);
for ii=3:10:2500
    M = imread([dir_,list(ii).name]); %raw unit8 images
    imagesc(M); 
    title(ii);
    drawnow; 
    pause();
end

%% find max position to locate LED

POS = [];
thre = 250;  %thrshold to find LED
for ii = 3:7:length(list)  %skipping through frames
    M = imread([dir_,list(ii).name]);
    [xx,yy] = find(M>thre);  %find LED
    POS = [POS; [mean(xx),mean(yy)]]; %center-of-mass of the LED
    ii
end

pix2mm = 1/31;  %pixel number to mm measured from the ruler image
POS = POS*pix2mm;
%% find PID signal
PID = [];
fileID = fopen([dir_,'PIDRead.txt'],'r');  %read value from PID
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
win = 30;  %window for averaging
del = 1*14;  %delay time after moving to that position
for ii = 1:25200 %3:7:length(list)-win-del-7  %subsample
    val = mean(A(ii+del:ii+del+win));
    PID = [PID; val]; %center-of-mass of the LED
    ii
end

%% plot
plate = imread('/tigress/LEIFER/Kevin/20190606/Data20190606_152844/Frame_000000.jpg');

figure()
ax1 = axes;
[x,y,z] = peaks;
imagesc(plate,'XData',[0 size(plate,2)*pix2mm],'YData',[0 size(plate,1)*pix2mm])
set(gca,'Ydir','reverse')
xlim([0 size(plate,2)*pix2mm])
ylim([0 size(plate,1)*pix2mm])
xlabel('X (mm)')
ylabel('Y (mm)');
view(2)
ax2 = axes;
scatter(POS(1:length(PID),2),POS(1:length(PID),1),30,PID,'filled')
set(gca,'Ydir','reverse')
xlim([0 size(plate,2)*pix2mm])
ylim([0 size(plate,1)*pix2mm])
%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
colormap(ax1,'gray')
colormap(ax2,'jet')

%% test interpolation
xx = POS(1:length(PID),1);
yy = POS(1:length(PID),2);
v = PID;
npos = find(isnan(xx)==1);
xx(npos) = [];
yy(npos) = [];
v(npos) = [];

[xq,yq] = meshgrid(min(xx):1:max(xx), min(yy):1:max(yy));
vq = griddata(xx,yy,v,xq,yq,'cubic');

mesh(xq,yq,vq)
hold on
plot3(xx,yy,v,'o')

%% Compare to SPG30
%SPG30
M1 = readtable('/tigress/LEIFER/Kevin/20190722/2019-07-22_13-27-05-SGP30_8487971.txt');
plot(M1.TVOC_SGP30_8487971)
SPG30 = M1.TVOC_SGP30_8487971;
%SPGC3
figure
M2 = readtable('/tigress/LEIFER/Kevin/20190722/2019-07-22_13-27-05-SGPC3_10997051.txt');
plot(M2.TVOC_SGPC3_10997051)
SPGC3 = M2.TVOC_SGPC3_10997051;
%%
%for raw signal
M = readtable('/tigress/LEIFER/Kevin/20190719/2019-07-19_14-22-45-raw.txt');
SGP30 = M.s_out_ETOH_SGP30_8487971;
SGPC3 = M.s_out_ETOH_SGPC3_10997336;
%plot
subplot(3,1,1:2)
plot(0.4*exp((0.5-SGP30)/512),'o')
hold on
plot(0.3*exp((0.5-SGPC3)/512),'ro')
subplot(3,1,3)
fac = length(SGPC3)/length(PID);
temp = [zeros(1,round(25))+0.0001 PID'];
plot(fac:fac:length(temp)*fac,log(temp))
%%
temp = SPG30(size(SPG30,1)-2*size(SPGC3,1):end);
nSPG30 = temp(1:2:2*length(SPGC3));
yyaxis right
plot(2:2:length(nSPG30)*2,nSPG30);
ylabel('SPG30')
yyaxis left
plot(2:2:length(SPGC3)*2,SPGC3)
ylabel('SPGC3 (ppb)')
xlabel('time (s)')
%%
%miniPID
PID = [];
fileID = fopen('/tigress/LEIFER/Kevin/20190703/PIDRead.txt','r');  %read value from PID
formatSpec = '%f';
A = fscanf(fileID,formatSpec);
for ii = 1:1:length(A)  %subsample
    val = A(ii);
    PID = [PID; val]; %center-of-mass of the LED
    ii
end

%%
%plot
hold on
yyaxis left
plot(PID);
yyaxis right
plot(SPG);
hold off
%%%%manually align...
subplot(211);plot(PID(100:end));xlim([0,length(PID(100:end))]);subplot(212);plot(SPG(50:end-60*9));xlim([0,length(SPG(50:end-60*9))])