%% Spatial distribution
%for angles
bin = 1;
filt = 10;
fract = [];
part = [1750,1500,1250,1000,750,500,250,0];
spart = cell(length(part)-1,1);  %spatial partition
cc = 1;
for c = 1:length(cand)
    
%     id = cand(c);
%     temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
%     temp1 = zeros(size(temp));
%     temp1(:,1) = smooth(temp(:,1),filt);
%     temp1(:,2) = smooth(temp(:,2),filt);
%     %temp1 = reshape(temp1,size(Tracks(id).Path,1),2);
%     temp2 = Tracks(id).Time;
        id = cand(c);
        temp2 = Tracks(id).Time;
%         plot(temp2); hold on
        overlap = find(temp2>0 & temp2<60*15);
        %%%only analyze overlaping time points
        if isempty(overlap)~=1
            
            
            %%%pre-processing
            temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
            temp1 = zeros(length(overlap),2);
            temp1(overlap,1) = smooth(temp(overlap,1),filt);
            temp1(overlap,2) = smooth(temp(overlap,2),filt);

    subs = temp1(1:bin:end,:);
    newtime = temp2(1:bin:end);
    vecs = [diff(subs); [0,0]];
    newtime = newtime(1:size(vecs,1));
    dists = zeros(1,size(subs,1));
    angs = zeros(1,size(vecs,1));    
    
    %%%select criteria
    p1 = subs(end,:); p2 = target;
    for dd = 1:length(dists)
        dists(dd) = distance(temp1(dd,:),target);
    end
    entp = find(dists<disth);
    %if distance(p1,p2) < disth  &&  distance(subs(1,:),p2) > disth  %&&   newtime(entp(1)) < 2000
    %if isempty(find(dists<disth)) ~= 1
    %if newtime(1)<60*20%  &&  newtime(1)> 60*5
    if 1==1
    
    %%%for distance
%     for dd = 1:length(dists)
%         dists(dd) = distance(subs(dd,:),target);
%     end
    infrac = zeros(1,size(subs,1));
    infrac(dists<disth) = 1;
    fract = [fract sum(infrac)/length(infrac)];
    trans = diff(infrac);
    tp = find(trans==1);

    
    %%%for angle
    for dd = 1:length(angs )
        angs(dd) = angles(vecs(dd,:),target-subs(dd,:));%ThetaInDegrees;%
    end
    
    %%%spatial dependency
    for pp = 1:length(part)-1
        pos = find(dists<part(pp) & dists>part(pp+1)); 
        spart{pp,cc} = angs(pos);%dists(pos);
    end
    cc = cc+1;
    
    end
    
        end
    
end

%% Spatial distribution
%for runs
bin = 1;
filt = 10;
fract = [];
part = [1750,1500,1000,500,0];
spart = cell(length(part)-1,1);  %spatial partition
cc = 1;
for c = 1:length(cand)
    
        id = cand(c);
        
        for pp = 1:length(part)-1
            
            temp = Tracks(id).Path;%smooth(Tracks(id).Path,6);
            temp2 = Tracks(id).Time;
            overlapt = find(temp2>0 & temp2<60*20);
            overlaps = find(temp(:,1)<part(pp) & temp(:,1)>part(pp+1)); 
            overlap = intersect(overlapt,overlaps);
            
        %%%only analyze overlaping time points
        if length(overlap)>=2
            
            %%%pre-processing
            temp1 = zeros(length(overlap),2);
            temp1(overlap,1) = smooth(temp(overlap,1),filt);
            temp1(overlap,2) = smooth(temp(overlap,2),filt);

            subs = temp1(1:bin:end,:);
            vecs = [diff(subs); [0,0]];
            dists = zeros(1,size(subs,1));
            angs = zeros(1,size(vecs,1));    
    
            %%%select criteria
            p1 = subs(end,:); p2 = target;
            for dd = 1:length(dists)
                dists(dd) = distance(temp1(dd,:),target);
            end
    %entp = find(dists<disth);
    %if distance(p1,p2) < disth  &&  distance(subs(1,:),p2) > disth  %&&   newtime(entp(1)) < 2000
    %if isempty(find(dists<disth)) ~= 1
    %if newtime(1)<60*20%  &&  newtime(1)> 60*5
    
            if 1==1
    
    %%%for distance
%     for dd = 1:length(dists)
%         dists(dd) = distance(subs(dd,:),target);
%     end

    infrac = zeros(1,size(subs,1));
    infrac(dists<disth) = 1;
    fract = [fract sum(infrac)/length(infrac)];
    trans = diff(infrac);
    tp = find(trans==1);

    
    %%%for angle
    for dd = 1:length(angs )
        angs(dd) = angles(vecs(dd,:),target-subs(dd,:));%ThetaInDegrees;%
    end
    
    %%%spatial dependency
%     for pp = 1:length(part)-1
%         pos = find(dists<part(pp) & dists>part(pp+1)); 
%         spart{pp,cc} = angs(pos);%dists(pos);
%     end
%     cc = cc+1;
    
    %%%run length
        [timestamps,runs] = def_turns(angs, 90, subs);
        spart{pp,cc} = runs;%dists(pos);
    
    cc = cc+1;
    
    end
    end
    
        end
    
end
%% spatial partition
figure;
for p = 1:size(spart,1)
    temp = [];
    for c = 1:size(spart,2)
        atr = spart{p,c};
        %if isempty(atr)~=1
        temp = [temp atr];
        %end
    end
    %figure; 
    %%%angle
%     subplot(1,size(spart,1),p); hist(temp); %title(num2str(part(p)*0.367)); 
    
    %%%run lengths
    [yy,xx] = hist(temp*0.0367,100);
    subplot(1,size(spart,1),p); semilogy(xx,yy,'r','LineWidth',3)
    xlim([0,15]); ylim([0,10000])
    
%     xlabel('angle'); ylabel('count')
%     xlabel('runs'); ylabel('count')
    set(gca,'Fontsize',20)
end



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


