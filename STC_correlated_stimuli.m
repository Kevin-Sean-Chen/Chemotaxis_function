%%%STC analysis
load('ADV.mat')

% tag = [];
% for i = 1:size(alltrigs,1); if sum(alltrigs(i,:))==0; tag = [tag i];end; end
% alltrigs(tag,:) = [];

sta = mean(alltrigs);
D = alltrigs-repmat(sta,size(alltrigs,1),1);
C = cov(D);
[uu,ss,vv] = svd(C);
PC1 = D*uu(:,1);
PC2 = D*uu(:,2);

%% C_prior
rX = [];
for ww = 1:size(Tracks,2)
win = 140;
acs = 140;
stim = Tracks(ww).LEDPower;
pos = randi([1,length(stim)],1,3);%find(temp(beh,:) == 1);

trigs = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(stim)
        trigs(ii,:) =stim(pos(ii)-win:pos(ii)+acs);
    end
end
rX = [rX; trigs];
end

sta_ = mean(rX);
D_ = rX - repmat(sta_, size(rX,1), 1);
Cp = cov(D_);

%%
dC = C-Cp;
[uu,ss,vv] = svd(dC);
PC1 = D*uu(:,1);
PC2 = D*uu(:,2);

kk = uu(:,1:2)\sta';

%%
Fs = 14;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 281;             % Length of signal
t = (0:L-1)*T;        % Time vector
Y = fft(smooth(mean(alltrigs),10));%(sin(0:281));%

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Self-defined correlated stimulus
freq = 0.2;
nt = 0.5;
x = 0:1/14:60*30;
y = sin(x*freq*(2*pi)) + nt*smooth(randn(1,length(x)),7,'moving')';
plot(x,y)

%%
y = y-mean(y);
y = y/std(y)*30;
y = y+30;
y(y>60) = 60;
y(y<0) = 0;
plot(x,y)

%%
fid = fopen('cor_stim_Kevin.txt','wt');
for ii = 1:length(y)
    fprintf(fid,'%g\t',y(ii));
    fprintf(fid,'  ');
end
fclose(fid)