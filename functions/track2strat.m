function [brw_index, wv_index] = track2strat(tracks, M)
%%% with simulated tracks, compute the emperical strategies (BRW and WV),
%%% calculated with Luo 2014 methods
%%% note that we ad-hoc define pirouettes and runs in this script

%%% parameters
dC_window = 10;  %time window for dC measurement for turns
turn_thre = 50;
tc_thre = 100;
min_x = 1; % remove little displacement segments like what we do for experiments(?
min_r = 1;

% initializing counts
run_up = 1;  %recording for runs
run_dn = 1;
turn_up = 1;  %recording for turns
turn_dn = 1;

for c = 1:length(tracks)
    
    %%% find runs
%     prs = find(abs(tracks(c).dth) > turn_thre);  % finding sharp turns
%     inter_pr = diff(prs);
%     pos = find(inter_pr < tc_thre);  % finding close consecutive turns
%     prs((pos)+1) = [];
    %%% testing burst detection
    prs = (abs(tracks(c).dth) > turn_thre);
    burst = conv(prs,ones(1,round(tc_thre)));%,'same');
    pos_init = find(diff(burst)>.5 & burst(1:end-1)==0);
    pos_fina = find(diff(burst)<-.5 & burst(2:end)==0);
    
    %%% serise of logic to align initial and final points of pr
    if length(pos_init)>2 && length(pos_fina)>2
        if length(pos_init)==length(pos_fina)
            if pos_init(1)>pos_fina(1)
                prs = [pos_init(2:end)'   pos_fina(1:end-1)'];
            else
                prs = [pos_init(1:end)'   pos_fina(1:end)'];
            end
        elseif length(pos_init)>length(pos_fina)
            prs = [];
            if pos_init(1)>pos_fina(1)
                prs = [pos_init(1:end-2)'   pos_fina(2:end)'];
            else
                prs = [pos_init(1:end-1)'   pos_fina(1:end)'];
            end
        elseif length(pos_init)<length(pos_fina)
            if pos_init(1)>pos_fina(1)
                prs = [pos_init(1:end)'   pos_fina(2:end)'];
            else
                prs = [];%[pos_init(1:end)'   pos_fina(1:end-1)'];
            end
%             if pos_init(1)>pos_fina(1)
%                 prs = [pos_init(1:end)'   pos_fina(2:end)'];
%             else
%                 prs = [];%[pos_init(1:end)'   pos_fina(1:end-1)'];  % ignore weird case for now...
%             end
        end
    else
        prs = [];
    end
        
%     try
%         prs = [pos_init(1:end)'   pos_fina(1:end)'];
%     catch
% %         prs = [pos_init'   pos_fina(1:end)'];
%         prs = [];
%     end
    %%%
    
    paths = tracks(c).xy';
    runs = [];
    if size(prs,1) > 2
        for ri = 1:length(prs)-1
%             runs = [runs; [prs(ri)  prs(ri+1)]];  % run are in between turns
              runs = [runs; [prs(ri,2)-dC_window*1  prs(ri+1,1)]];   % the last pr end to next pr on
        end
    end
    
    %%% computing runs (biased-random walk)
    if isempty(runs)==0
    for rr = 1:size(runs,1)
        path_i = paths(runs(rr,1):runs(rr,2),:);  %segment of running
        displace = mean((path_i(:,1)-mean(path_i(:,1))).^2 + (path_i(:,2)-mean(path_i(:,2))).^2); %pixel displacement
        if displace > min_x^2 && length(runs(rr,1):runs(rr,2)) > min_r  % only compute for reasonable runs
            dC = -(M(floor(path_i(1,2)),floor(path_i(1,1))) - M(floor(path_i(end,2)),floor(path_i(end,1))));
            if dC>=0  %up gradient  %%%%%%%%% hack for now to invert concentration~~~
                run_up = [run_up    length(path_i)];%sqrt((path_i(1,1)-path_i(end,1))^2+(path_i(1,2)-path_i(end,2))^2)];%*dC];  %record run length
            elseif dC<0
                run_dn = [run_dn    length(path_i)];%sqrt((path_i(1,1)-path_i(end,1))^2+(path_i(1,2)-path_i(end,2))^2)];% length(path_i)];%*abs(dC)]; 
            end      
        end
    end
    end
    

    %%% find turns
%     dC_window = 10;%10*5/14;
%     prs = find(abs(tracks(c).dth) > 100);  % finding sharp turns
%     inter_pr = diff(prs);
%     pos = find(inter_pr < 100);  % finding close consecutive turns, with longer interval here.
%     prs((pos)+1) = [];
    
    %%% computing turns (weathervaning)
    if isempty(prs)==0
    for tt = 1:size(prs,2)-1
%         c_i = max((prs(tt)-round(dC_window)*1), 1);  %initiation time 
%         c_f = min((prs(tt+1)+round(dC_window)*1), length(paths));  %exiting time...  note that in sumulation we don't have an explicit piroutte state!!!
        c_i = max((prs(tt,1)-round(dC_window)*0), 1);  %initiation time 
        c_f = min((prs(tt,2)-round(dC_window)*1), length(paths)); 
        
        dC = -(M(floor(paths(c_i,2)),floor(paths(c_i,1))) - M(floor(paths(c_f,2)),floor(paths(c_f,1))));
        if dC>=0  %up gradient
            turn_up = turn_up + 1;%dC;  %record run length
        elseif dC<0
            turn_dn = turn_dn + 1;%abs(dC); 
        end   
    end
    end
    
end

%%% analysis, Luo 2014 methods
brw_index = (median(run_up)-median(run_dn)) / (median(run_up)+median(run_dn));
wv_index = ((turn_up)-(turn_dn)) / ((turn_up)+(turn_dn));

end