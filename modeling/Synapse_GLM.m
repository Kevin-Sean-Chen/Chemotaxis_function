%%%Synapse_GLM
%two synaptically coupled spliking neurons driven by noisy input
%coupled GLM fitting for this simple circuit
%% stimuli
T = 8000;
dt = 0.1;
time = dt:dt:T;
ll = length(time);
Ast = 30;
stim = randn(1,ll)*Ast;

%% spiking neurons
%Izhikevich parameters
a = 0.02;
b = 0.2;
c = -55;
d = 4.;
%synaptic paramters
U = 0.5;
tauD = 0.2;
tauF = 1.5;
Asyn = 400;

%initialize dynamical variables
vs = zeros(2,ll);  %neural
ws = zeros(2,ll);
us = zeros(1,ll);  %synapse
xs = zeros(1,ll);
us(1) = U;
xs(1) = 1;
Gs = zeros(1,ll);
spk1 = zeros(1,ll);  %spikes trains
spk2 = zeros(1,ll);

%noise sourse
v1n = 0;
v2n = 0;
gn = 0;

%looping
for t = 1:ll-1
    %%%pre-synaptic neuron
    vs(1,t+1) = vs(1,t) + dt*(0.04*vs(1,t)^2 + 5*vs(1,t) + 140 - ws(1,t) + stim(t)) + dt^0.5*randn*v1n;
    ws(1,t+1) = ws(1,t) + dt*(a*(b*vs(1,t) - ws(1,t)));
    if vs(1,t+1)>=30
        vs(1,t+1) = c;
        ws(1,t+1) = ws(1,t) + d;
        spk1(t+1) = 1;
    end
    
    %%%synapse
    xs(t+1) = xs(t) + dt*((1-xs(t))/tauD + us(t)*xs(t)*spk1(t)*Asyn);
    us(t+1) = us(t) + dt*((U-us(t))/tauF + U*(1-us(t))*spk1(t)*Asyn);
    Gs(t+1) = max(xs(t)*us(t) + gn*randn,0);
    
    %%%post-synaptic neuron
    vs(2,t+1) = vs(2,t) + dt*(0.04*vs(2,t)^2 + 5*vs(2,t) + 140 - ws(2,t) + Gs(t)) + dt^0.5*randn*v2n;
    ws(2,t+1) = ws(2,t) + dt*(a*(b*vs(2,t) - ws(2,t)));
    if vs(2,t+1)>=30
        vs(2,t+1) = c;
        ws(2,t+1) = ws(2,t) + d;
        spk2(t+1) = 1;
    end
end

%% STA!
win = 150;
X = zeros(ll-win,win);
% for i = 1:ll-win
%     X(i,:) = stim(i:i+win-1);
% end
acs = 30;
X = zeros(ll-win,win+acs);
for i = acs+1:ll-win
    X(i-acs,:) = stim(i-acs:i+win-1);
end
STA1 = (X'*X)\X'*vs(1,win+1:end)';
STA2 = (X'*X)\X'*vs(2,win+1:end)';

%% LL optimization (minimize negative log-likelihood)
% fun = @(x) -vs(1,win+1:end)*log(X*x') + sum(X*x');
% x0 = ones(1,win+acs);
% options = optimset('Maxiter',1000);
% x = fminsearch(fun,x0,options);

%% Coupling!
coup = 60;
Xc = zeros(ll-coup,coup);
for i = 1:ll-coup
    Xc(i,:) = vs(1,i:i+coup-1);
end
coupK = (Xc'*Xc)\Xc'*vs(2,coup+1:end)';

%% spike history
spkh = 50;
xh = zeros(ll-spkh,spkh);
for i = 1:ll-spkh
    xh(i,:) = vs(2,i:i+spkh-1);
end
spk_hist = (xh'*xh)\xh'*vs(2,spkh+1:end)';

%% %%% plots %%%
%% activity
% figure;
% plot(vs(1,:))
% hold on
% plot(vs(2,:))
% plot(stim*0.1,'k')
% %% STA
% figure;
% plot(STA1)
% hold on
% plot(STA2)
% %% Coupling
% figure;
% plot((-coup)*dt:dt:-dt,coupK)
% %% History
% figure;
% plot(spk_hist)

%% scanning
% uus = 0.1:0.1:0.9;
% cc = hsv(length(uus));
% for kk = 1:length(uus); U = uus(kk);  Synapse_GLM; plot((-coup)*dt:dt:-dt,coupK,'color',cc(kk,:)); hold on; end


