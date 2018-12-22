%Check_frames
%cd 'Z:\Ali\data\group2\Data20180309_160703_Naive1'
cd 'Z:\Ali\data\group2-new_lights\Data20180309_160703_Naive1'
listing = dir('Frame_*.jpg');

nf = size(listing,1);
vals = [];
for i = 1:nf
    temp1 = listing(i,1).name;
    temp2 = extractBetween(temp1,'Frame_','.jpg');
    nn = str2num(temp2{1});
    vals = [vals nn];
end

plot(vals)