%% Common traces
%
bin = 1;
filt = 6;
fract = [];
cendtr = {};  %common ending traces
cc = 1;
for c = 1:length(cand)
    
    id = cand(c);
    temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
    temp1 = zeros(size(temp));
    temp1(:,1) = smooth(temp(:,1),filt);
    temp1(:,2) = smooth(temp(:,2),filt);
    %temp1 = reshape(temp1,size(Tracks(id).Path,1),2);
    temp2 = Tracks(id).Time;
    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [diff(subs); [0,0]];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    p1 = subs(end,:); p2 = target;
%     for dd = 1:length(dists)
%         dists(dd) = distance(temp1(dd,:),target);
%     end
    %if distance(p1,p2) < disth  &&  distance(subs(1,:),p2) > disth
    if isempty(find(dists<disth)) ~= 1
    
    %%%for distance
    for dd = 1:length(dists)
        dists(dd) = distance(subs(dd,:),target);
    end
    infrac = zeros(1,size(subs,1));
    infrac(dists<disth) = 1;
    fract = [fract sum(infrac)/length(infrac)];
    trans = diff(infrac);
    tp = find(trans==1);

    
    %%%for angle
    for dd = 1:length(angs)
        angs(dd) = angles(vecs(dd,:),target-subs(dd,:));%ThetaInDegrees;%
    end
    
    
    %%%condition for entering zone
    if isempty(tp)~=1
        %tp = tp(end-max(0,length(tp)-3):end);
        length(tp);
        %tp = tp(end);
        for tps = 1:length(tp)
            cendtr{cc} = angs(1:tp(tps));
            cc = cc+1;
        end
    end
    
    end
    
end

%% collapse common traces
bk = 1000;
altr = [];
for i = 1:length(cendtr)
    if length(cendtr{i}) <= bk
        continue
    else
        altr = [altr;cendtr{i}(end-bk:end)];
    end
end

errorbar(-fliplr([1:bk+1]./14),mean(altr),std(altr))
hold on
plot(-fliplr([1:bk+1]./14),mean(altr),'ro')
xlabel('time before reaching')
ylabel('angle')

% %functions in the end~
% function dd = distance(p1,p2)%p2=target
%     dd = sqrt((p1(1)-p2(1))^2 + (p1(2)-p2(2))^2);
% end
% function aa = angles(v,u)%u=target
%     CosTheta = dot(u,v)/(norm(u)*norm(v));
%     ThetaInDegrees = acosd(CosTheta);
%     aa = ThetaInDegrees;
% end
% function nt = normtime(time)
%     nt = time-time(1);
% end


