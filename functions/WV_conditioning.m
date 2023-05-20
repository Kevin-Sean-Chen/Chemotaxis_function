%% WV_conditioning
%%%%
% To analyze weathervaning using conditions on the odor concentration
% This is more comparible to Iino 2009 analysis
% Tracks and cand are structure and vectors available
%%%%

x_condition = 0;  % try a line in space first
filt = 14*2.;
fixlen = 1;  %we only compute a 2d vector when it moves forward this much
n_smooth = 200;
[fx,fy] = gradient((conv2(M,ones(n_smooth,n_smooth),'same')/n_smooth^2),1);  %prepare 2D gradient field
all_grad_ = [];
all_angs_ = [];
all_amps_ = [];  % condition on gradient amplitude

figure()
% imagesc(M); hold on
for nn = 1:length(Tracks)
    if isempty(find(cand==nn))==0

vec_i = [];    %2d vectors with adaptive sampling rate
grad_i = [];   %at this location, what was the gradient direction
sped_i = [];

speeds = Tracks(nn).SmoothSpeed;
paths = Tracks(nn).Path;
angspeed_i = Tracks(nn).AngSpeed;
runs = Tracks(nn).Runs;

for rr = 1:size(runs,1)
        path_i = paths(runs(rr,1):runs(rr,2),:); %paths; %
        speed_i = speeds(runs(rr,1):runs(rr,2)); 
        x_smooth = smooth(path_i(:,1), filt,'sgolay',poly_degree);
        y_smooth = smooth(path_i(:,2), filt,'sgolay',poly_degree);
        path_i = [x_smooth   y_smooth];  
        
ii = 1;
pos_temp = path_i(1,:);  %initial location
% if M(floor(path_i(1,2)),floor(path_i(1,1))) < M(floor(path_i(end,2)),floor(path_i(end,1))) %&& path_i(end,1)<x_condition
while ii<length(speed_i)
    
%     delta_t = min(floor(fixlen/speed_i(ii)), floor(length(speed_i)/2));  %discrete time window to update, adpated to the velocity
    delta_t = min(14*fixlen, floor(length(speed_i)/2));%
    vec_i = [vec_i; path_i(ii,:)-pos_temp];  %compute dispacement vector
    grad_i = [grad_i; [fx(floor(path_i(ii,2)),floor(path_i(ii,1))), fy(floor(path_i(ii,2)),floor(path_i(ii,1)))] ];  %gradident direction at this location
%     grad_i = [grad_i; ([0,2500] - path_i(ii,:)) ];  % testing global bearing
    pos_temp = path_i(ii,:);  %update postion
    sped_i = [sped_i; speed_i(ii)];
    ii = ii + delta_t+1;
    
% end
end

angs = zeros(1, length(vec_i)-1);
grad = angs*1;
amps = angs*1;
for pp = 1:length(vec_i)-1
    angs(pp) = -angles(vec_i(pp+1,:)/norm(vec_i(pp+1,:)),vec_i(pp,:)/norm(vec_i(pp,:))) / (norm(vec_i(pp,:))*pix2mm); %/fixlen;% / ((sped_i(pp)+sped_i(pp+1))/2);%   %
    grad(pp) = angles(vec_i(pp,:)/norm(vec_i(pp,:)), grad_i(pp,:)/norm(grad_i(pp,:)));% *norm(grad_i(pp,:)) * pix2mm;  % bearing
    
%     grad(pp) = sin((angles(vec_i(pp,:)/norm(vec_i(pp,:)), grad_i(pp,:)/norm(grad_i(pp,:)))) / pi) * norm(grad_i(pp,:));  % measure C^\perp
%     grad(pp) = dot(vec_i(pp,:), grad_i(pp,:));%
%     grad(pp) = angles(vec_i(pp,:)/norm(vec_i(pp,:)), [1,0]);% *norm(grad_i(pp,:));  % mock x-axis flow direction
    amps(pp) = norm(grad_i(pp,:)); %sin(angles(vec_i(pp,:)/norm(vec_i(pp,:)), grad_i(pp,:)/norm(grad_i(pp,:)))/pi)*   
end

all_grad_ = [all_grad_  grad];
all_angs_ = [all_angs_  angs];
all_amps_ = [all_amps_  amps];
nn

end

% plot(path_i(:,1),path_i(:,2));  hold on
% subplot(131); plot(grad); hold on; plot(angs);plot(angspeed_i); hold off;
% % subplot(121); 
% subplot(132); plot(path_i(:,1),path_i(:,2))
% subplot(133); plot(speed_i)
% pause();

    end
end

%%
% pass variables for conditioning
all_grad = all_grad_;
all_angs = all_angs_;
all_amps = all_amps_;

%% conditions
% pos = find(all_grad==0 | all_grad==0);
% all_grad(pos) = [];
% all_angs(pos) = [];

pos = find(abs(all_amps)<0.04);%nanstd(all_amps)+nanmean(all_amps));% | abs(all_angs)>1000);  %find(abs(all_angs)>200); %
all_grad(pos) = NaN;
all_angs(pos) = NaN;

%% WV!
nbs = linspace(-150,150,15);%
thr = nanstd(all_grad)*2
figure()
% pos_0 = find((all_grad)<-thr); %
% pos_60 = find(abs(all_grad)<thr);
% pos_120 = find((all_grad)>thr);

pos_0 = find((all_grad)>-135 & all_grad<-45);
pos_60 = find(abs(all_grad)<45 | abs(all_grad)>135);
pos_120 = find((all_grad)>45 & all_grad<135);
% pos_test = find(abs(all_grad)>135);

% pos_0 = find((all_grad)>-110 & all_grad<-70);
% pos_60 = find(abs(all_grad)<20 | abs(all_grad)>160);
% pos_120 = find((all_grad)>70 & all_grad<110);

% pos_0 = find((all_grad)>-120 & all_grad<-60);
% pos_60 = find(abs(all_grad)<30);
% pos_120 = find((all_grad)>60 & all_grad<120);

H1 = histogram(all_angs(pos_0), nbs, 'Normalization', 'pdf'); hold on
H2 = histogram(all_angs(pos_60), nbs, 'Normalization', 'pdf'); hold on
H3 = histogram(all_angs(pos_120), nbs, 'Normalization','pdf'); hold on
% H4 = histogram(all_angs(pos_test), nbs, 'Normalization','pdf');
% close(fig);

figure()
aa = H1.Values;  bb = H1.BinEdges;
bb = (bb(2:end) + bb(1:end-1))/2;
plot(bb,aa/sum(aa) / 1); hold on %max(aa/sum(aa))
med = mean((aa/sum(aa) / max(aa/sum(aa))).*bb);
plot([med,med],[0,1])
y = skewness(all_angs(pos_0))
hold on
aa = H2.Values;  bb = H2.BinEdges;
bb = (bb(2:end) + bb(1:end-1))/2;
plot(bb,aa/sum(aa) / 1); hold on
med = mean((aa/sum(aa) / max(aa/sum(aa))).*bb);
plot([med,med],[0,1])
y = skewness(all_angs(pos_60))
hold on
aa = H3.Values;  bb = H3.BinEdges;
bb = (bb(2:end) + bb(1:end-1))/2;
plot(bb,aa/sum(aa) / 1); hold on
med = mean((aa/sum(aa) / max(aa/sum(aa))).*bb);
plot([med,med],[0,1])
y = skewness(all_angs(pos_120))

% aa = H4.Values;  bb = H4.BinEdges;
% bb = (bb(2:end) + bb(1:end-1))/2;
% plot(bb,aa/sum(aa) / 1); hold on
% med = mean((aa/sum(aa) / max(aa/sum(aa))).*bb);
% plot([med,med],[0,1])
% y = skewness(all_angs(pos_test))

xlabel('curvature (degrees/mm)'); ylabel('density')
set(gca,'Fontsize',20);  set(gcf,'color','w'); ylim([0,0.25])
%% WV tuning analysis (gradient vs. curvature)

pos = find(abs(all_angs)>180); %find(all_amps<0.05);
all_grad(pos) = NaN;
all_angs(pos) = NaN;

%%
figure
bins = 7;
downsamp = 5;
H = histogram(all_grad, bins);
bb = H.BinEdges;
% bb(end) = bb(end) + 0*H.BinWidth;
muv = zeros(1,bins);
epsv = zeros(2,bins);

for bi = 2:bins+1%+2
    pos = find(all_grad>bb(bi-1) & all_grad<bb(bi));
    temp0 = all_angs(pos);
    temp = nanmean(reshape([temp0(:); nan(mod(-numel(temp0),downsamp),1)],downsamp,[]));
    muv(bi-1) = nanmean(temp);%median(temp); %
    epsv(:,bi-1) = nanstd(temp)/sqrt(length(pos));% quantile(temp,[.25 .75],2);  %
end
bb = bb+H.BinWidth*.5;

%%

figure()
% plot(all_grad, all_angs,'.','color', [.7 .7 .7])
hold on; errorbar(bb(1:end-1), muv, 0+epsv(1,:) ,0+epsv(2,:),  'k', 'Linewidth',1) %+mean(diff(bb(1:end-1)))/1
% hold on; errorbar(bb(1:end-1), muv, muv-epsv(1,:) ,muv + epsv(2,:),  'k', 'Linewidth',5) %+mean(diff(bb(1:end-1)))/1
hold on; plot(bb(1:end-1), muv, 'k-o', 'Linewidth',1); 
%plot([min(all_delc),max(all_delc)],[0,0],'--','color', [.5 .5 .5])
xlabel('bearing to local gradient'); ylabel('curvature (degrees/mm)'); set(gcf,'color','w'); set(gca,'FontSize',20)
% ylim([-100,100])


%%
figure;
bins = 5;
num_points = floor(length(all_grad)/bins);  %number of points in one bin, now that we fix this value
[val,poss] = sort(all_grad);
bb = zeros(1,bins+1);
muv = zeros(1,bins);
epsv = zeros(2,bins);

for bi = 2:bins
%     pos = find(all_delc>bb(bi-1) & all_delc<bb(bi));
    pos = poss((bi-1)*num_points+1:(bi)*num_points);
    bb(bi) = nanmedian(all_grad(pos));
    temp0 = all_angs(pos);
    temp = nanmean(reshape([temp0(:); nan(mod(-numel(temp0),downsamp),1)],downsamp,[]));
    muv(bi-1) = nanmean(temp);%median(temp); %
    epsv(:,bi-1) = nanstd(temp)/sqrt(length(pos));% quantile(temp,[.25 .75],2);  %
end
