function [brw_index, wv_index] = track2strat(tracks, M)
%%% with simulated tracks, compute the emperical strategies (BRW and WV),
%%% calculated with Luo 2014 methods

%%% parameters
dC_window = 10;  %time window for dC measurement for turns
turn_thre = 100;
tc_thre = 2;
min_x = 10; % remove little displacement segments like what we do for experiments(?

% initializing counts
run_up = 1;  %recording for runs
run_dn = 1;
turn_up = 1;  %recording for turns
turn_dn = 1;

for c = 1:length(tracks)
    
    %%% find runs
    prs = find(abs(tracks(c).dth) > turn_thre);  % finding sharp turns
    inter_pr = diff(prs);
    pos = find(inter_pr < tc_thre);  % finding close consecutive turns
    prs((pos)+1) = [];
    paths = tracks(c).xy';
    runs = [];
    if length(prs) > 1
        for ri = 1:length(prs)-1
            runs = [runs; [prs(ri)  prs(ri+1)]];  % run are in between turns
        end
    end
    
    %%% computing runs (biased-random walk)
    if isempty(runs)==0
    for rr = 1:size(runs,1)
        path_i = paths(runs(rr,1):runs(rr,2),:);  %segment of running
        displace = mean((path_i(:,1)-mean(path_i(:,1))).^2 + (path_i(:,2)-mean(path_i(:,2))).^2); %pixel displacement
        if displace > min_x^2   % only compute for reasonable runs
            dC = -(M(floor(path_i(1,2)),floor(path_i(1,1))) - M(floor(path_i(end,2)),floor(path_i(end,1))));
            if dC>=0  %up gradient  %%%%%%%%% hack for now to invert concentration~~~
                run_up = [run_up length(path_i)];%*dC];  %record run length
            elseif dC<0
                run_dn = [run_dn length(path_i)];%*abs(dC)]; 
            end      
        end
    end
    end
    

    %%% find turns
    dC_window = 10;%10*5/14;
    prs = find(abs(tracks(c).dth) > 100);  % finding sharp turns
    inter_pr = diff(prs);
    pos = find(inter_pr < 100);  % finding close consecutive turns
    prs((pos)+1) = [];
    
    %%% computing turns (weathervaning)
    if isempty(prs)==0
    for tt = 1:size(prs,2)-1
        c_i = max((prs(tt)-round(dC_window)*0), 1);  %initiation time 
        c_f = min((prs(tt+1)+round(dC_window)*0), length(paths));  %exiting time...  note that in sumulation we don't have an explicit piroutte state!!!
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
brw_index = (mean(run_up)-mean(run_dn)) / (mean(run_up)+mean(run_dn));
wv_index = ((turn_up)-(turn_dn)) / ((turn_up)+(turn_dn));

end