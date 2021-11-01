%Chemotaxis_Theory

%% Convolution effects
t = 0:1:300;
K = a*exp(-a*t).*((a*t).^5/factorial(5) - (a*t).^7/factorial(7));

x = randn(1,5000);
y1 = smooth(x,find(K == max(K)));
y2 = smooth(x,find(K == min(K)));

figure;
plot(conv(K,y1))
hold on
plot(conv(K,y2))

%% R&T strategy
