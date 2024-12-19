%ARprocess_test
%this is a test for modeling head angles with different history structures
%% AR(1)
dt = 0.1;
T = 1000;
time = 0:dt:T;
x = zeros(1,length(time));
x(1) = randn(1);

phi1 = 2.5;
sig1 = 0.5;

for tt = 1:length(x)-1
    x(tt+1) = phi1*x(tt)*dt + sig1*randn(1)*sqrt(dt);  %adding noise
end

figure()
plot(time,x)

figure()
nlags = 100;
xx = -nlags:nlags;
xcsamp = xcov(x-mean(x),nlags, 'unbiased');
plot(xx,xcsamp, '.-');


%% AR(2)
dx = zeros(1,length(time));
dx(1) = randn(1);
phi2 = 1.5;
sig2 = 0.5;
for tt = 1:length(x)-1
    x(tt+1) = x(tt) + dx(tt);  %adding noise to change
    dx(tt+1) = phi2*dx(tt)*dt + sig2*randn(1)*sqrt(dt);
end

figure()
plot(time,x)

figure()
nlags = 100;
xx = -nlags:nlags;
xcsamp = xcov(x-mean(x),nlags, 'unbiased');
plot(xx,xcsamp, '.-');
