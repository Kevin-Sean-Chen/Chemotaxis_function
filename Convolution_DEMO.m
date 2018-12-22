%Convolution_DEMO
X = Tracks(1).LEDPower;
Ker = mean(alltrigs);
Ker = Ker-mean(Ker);
% X = conv(X,Ker);
% X = X(length(Ker):end);
%X = Ker;

Fs = 14;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(X);             % Length of signal
t = (0:L-1)*T;        % Time vector
Y = fft(X);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
P1 = P1/sum(P1);
plot(f,P1,'LineWidth',3) 
xlabel('f (Hz)')
ylabel('PSD')
xlim([0,1])
set(gca,'Fontsize',30)