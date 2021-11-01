%%%RL_test
%load('/tigress/LEIFER/Ali/data/group3-cam_adjust/Data20180330_141117_Ap4/analysis/Path.mat')
%load('/Data20180330_141117_Ap4/Path.mat')
%load('/tigress/LEIFER/Ali/data/group3-cam_adjust/Data20180330_131220_N4/analysis/Path.mat')

%% loading data


%% Embedding
temp = values;
temp2 = temp;
alltr = [];
temp_p = [];
temp_f = [];
for ii = 1:length(temp)
    temp2{ii} = temp{ii}(:,1)-temp{ii}(1,1);
    %alltr = [alltr temp2{ii}'];
    temp_p = [temp_p temp2{ii}(1:end-1)'];
    temp_f = [temp_f temp2{ii}(2:end)'];
end

%% Transition matrix
Nbs = 20;
temp_p2 = downsample(temp_p,10);
temp_f2 = downsample(temp_f,10);
[Joint,c] = hist3([temp_p2; temp_f2]','Nbins',[Nbs,Nbs]);
[margin_p,centers] = hist(temp_p2,c{1});
[margin_f,centers] = hist(temp_f2,c{2});
Joint(Joint==0) = 1; margin_p(margin_p==0) = 1; margin_f(margin_f==0) = 1;
Joint = Joint/sum(sum(Joint));
margin_p = margin_p/sum(margin_p);
margin_f = margin_f/sum(margin_f);
Trans = Joint./margin_p;
Trans = Trans./sum(Trans);

%% Minimization
lbd = 10;
%passive dynamic P(s'|s)
xx = -Nbs/2:Nbs/2;
D = .5;
diffu = exp(-(xx-0).^2/(2*D^2));
pss = rand(Nbs,Nbs) + (1-diffu(2:end)*2);
pss = pss./sum(pss,1);
f = @(x) -sum(log(x)'*margin_f') + sum(log(pss*x)'*margin_p') + lbd*sum(diff(log(x))'.^2);
x0 = rand(size(pss,1),1)*5;
options = optimset('MaxFunEvals',10000)
options = optimset('MaxIter',10000)
x = fminsearch(f,x0,options);

plot(centers,-log(x))
