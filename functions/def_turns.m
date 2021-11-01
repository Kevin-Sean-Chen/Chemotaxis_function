function [timestamps,runs] = def_turns(angles,threshold,subs)

da = angles;%diff(angles);
[val,pos] = find(da>threshold);% & da<170);
timestamps = pos;

runs = zeros(1,length(pos)-1);
for dd = 1:length(runs)
    runs(dd) = distance(subs(timestamps(dd),:),subs(timestamps(dd+1),:));
end


end