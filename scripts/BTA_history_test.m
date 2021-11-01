%BTA_history_test
bb = 8%1:size(Tracks(1).Behaviors,1)
alltrigs = [];
allhist = [];

for ww = 1:size(Tracks,2)
num = ww;
beh = bb; % behavioral state
win = 200;
acs = 200;
temp = Tracks(ww).Behaviors;
stim = Tracks(ww).LEDPower;
pos = find(temp(beh,:) == 1);
%%%if isempty(pos)~=1; pos = randi(length(temp(beh,:)),length(pos)); end
%pos = find(sum(allB{num}) >= 1);
trigs = zeros(length(pos),win+acs+1);
hists = zeros(length(pos),win+acs+1);
for ii = 1:length(pos)
    if pos(ii)-win>0 && pos(ii)+acs<length(temp(beh,:))
        trigs(ii,:) =stim(pos(ii)-win:pos(ii)+acs);
        hists(ii,:) =temp(beh,pos(ii)-win:pos(ii)+acs);
    end
end
alltrigs = [alltrigs; trigs];
allhist = [allhist; hists];
end
% %figure; plot([-acs:win]*(1/14),mean(alltrigs))
% set(gca, 'XTickLabel', [],'XTick',[])
% subplot(size(Tracks(1).Behaviors,1),1,bb); plot([-acs:win]*(1/14),mean(alltrigs))
% set(gca,'Fontsize',20)