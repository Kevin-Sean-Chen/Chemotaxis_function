%head_analysis
load('/projects/LEIFER/Sandeep/Learning/GLM_KSC_SK/tracks_of_interest_AML67_20200902_165911.mat')
sample = tracks_of_interest(10);

%% compute head anlge
lt = length(sample.Frames);  %time vector
head = [2,4];  %define head segment
neck = [4,7];  %define neck segment
d_theta = zeros(1,lt);
time = [1:lt]*1/30;

for tt = 1:lt
    cls = sample.Centerlines(:,:,tt);
    vec_head = cls(head(1),:) - cls(head(2),:);
    vec_body = cls(neck(1),:) - cls(neck(2),:);
    d_theta(tt) = angles(vec_head, vec_body);  %measure head angles
end

%% plotting
figure;
subplot(3,2,1)
plot(time,d_theta);ylabel('head angle'); xlim([0,max(time)]); set(gca,'FontSize',16)
subplot(3,2,3)
plot(time,sample.EllipseRatio);ylabel('ellipseratio'); xlim([0,max(time)]); set(gca,'FontSize',16)
subplot(3,2,5)
imagesc(sample.BehavioralAnnotation); ylabel('behavior')
set(gca, 'XTick', [0:0.2:1]*lt, 'XTickLabel', [0:0.2:1]*max(time)) % 10 ticks 
ylabel('time (s)'); set(gca,'FontSize',16)
subplot(3,2,[2,4,6])
pos = find(sample.BehavioralAnnotation ~= 1);
plot(sample.Path(:,1),sample.Path(:,2)); hold on
plot(sample.Path(pos,1),sample.Path(pos,2),'ro')
set(gcf,'color','w'); set(gca,'FontSize',16)


%% simple spectral check for head swings
Fs = 30;  %sample rate
L = length(d_theta);
Y = fft(d_theta);
f = Fs*(0:(L/2))/L;
plot(f,Y) 
title('head swing spectrum')
xlabel('f (Hz)')
ylabel('|P1(f)|')
set(gcf,'color','w'); set(gca,'FontSize',20)

%% function for angle calculation
function aa = angles(v,u)%u=target
    %%%with angle sign
    v_3d = [v, 0];
    u_3d = [u, 0];
    c = cross(v_3d, u_3d);
    
    % calculate degrees
    if c(3) < 0
        aa = -atan2(norm(c),dot(v,u))*(180/pi);
    else
        aa = atan2(norm(c),dot(v,u))*(180/pi);
    end
end