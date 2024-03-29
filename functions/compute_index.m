function [ci_, brw_index, wv_index] = compute_index(Tracks, M, T);

%%% parameters
dC_window = 14*10;  %time window for dC measurement for turns
time_wind = 60*T;  %first few minutes
min_t = 60*2;  % minumum time window
min_x = 100;  % min displacement pixels
conc_thre = .6*max(max(M));  % high concentration initial to be excluded (60%)

% initializing counts
run_up = 0;  %recording for runs
run_dn = 0;
turn_up = 0;  %recording for turns
turn_dn = 0;
ci_low = 0;  %recording for chemotaxis performance
ci_high = 0;
ci1 = 1;  ci2 = 1; 

for c = 1:length(Tracks)
    %%% load observations
    prs = Tracks(c).Pirouettes;
    paths = Tracks(c).Path;
    runs = Tracks(c).Runs;
    times = Tracks(c).Time;
    
    %%%%% conditioning tracks (in space or in time)
    displace = mean((Tracks(c).Path(:,1)-mean(Tracks(c).Path(:,1))).^2 + (Tracks(c).Path(:,2)-mean(Tracks(c).Path(:,2))).^2); %pixel displacement
    if max(times)<=time_wind && Tracks(c).Time(end)-Tracks(c).Time(1) > min_t  && displace>min_x^2  &&  M(floor(Tracks(c).Path(1,2)), floor(Tracks(c).Path(1,1))) < conc_thre
    
    %%% computing runs (biased-random walk)
    if isempty(runs)==0
    for rr = 1:size(runs,1)
        path_i = paths(runs(rr,1):runs(rr,2),:);  %segment of running
%         dC = Fcon(path_i(1,1) , path_i(1,2)) - Fcon(path_i(end,1) , path_i(end,2));  %delta C
        dC = -(M(floor(path_i(1,2)),floor(path_i(1,1))) - M(floor(path_i(end,2)),floor(path_i(end,1))));
        if dC>=0  %up gradient  %%%%%%%%% hack for now to invert concentration~~~
            run_up = [run_up length(path_i)];%*dC];  %record run length
        elseif dC<0
            run_dn = [run_dn length(path_i)];%*abs(dC)]; 
        end      
    end
    end
    
    %%% computing turns (weathervaning)
    if isempty(prs)==0
    for tt = 1:size(prs,2)-1 %size(prs,1) %
        c_i = max((prs(tt,1)-1*dC_window), 1);  %initiation time 
        c_f = min((prs(tt,2)+1*dC_window), length(paths));  %exiting time
%         dC = Fcon(paths(c_i,1) , paths(c_i,2)) - Fcon(paths(c_f,1) , paths(c_f,2));  %delta C
        dC = -(M(floor(paths(c_i,2)),floor(paths(c_i,1))) - M(floor(paths(c_f,2)),floor(paths(c_f,1))));
        if dC>=0  %up gradient
            turn_up = turn_up + 1;%dC;  %record run length
        elseif dC<0
            turn_dn = turn_dn + 1;%abs(dC); 
        end   
    end
    end
    
    %%% quantify chemotaxis index
%     dC = Fcon(paths(1,1) , paths(1,2)) - Fcon(paths(end,1) , paths(end,2)); 
    c0 = M(floor(paths(1,2)),floor(paths(1,1)));
    dC = -(M(floor(paths(1,2)),floor(paths(1,1))) - M(floor(paths(end,2)),floor(paths(end,1))));
    if dC>=0
        ci_high = ci_high+1;%dC;  %recording for chemotaxis performance
% % % % % % %         ci_high = ci_high + dC;
%             ci_high = ci_high + dC/(100-c0);
%             ci1 = ci1 + 1;
    elseif dC<0
        ci_low = ci_low+1;%abs(dC);
% % % % % % %         ci_low = ci_low + dC;
%             ci_low = ci_low + dC/(c0-10);
%             ci2 = ci2 + 1;
    end
    
    end
    
end

% %% analysis
brw_index = (mean(run_up)-mean(run_dn)) / (mean(run_up)+mean(run_dn));
wv_index = ((turn_up)-(turn_dn)) / ((turn_up)+(turn_dn));
ci_ = (ci_high - ci_low) / (ci_high+ci_low);
% % % % % % % ci_ = (ci_high + ci_low) / (ci_high+abs(ci_low));
% ci_ = (ci_high/ci1 + ci_low/ci2);

% disp(folder_names);
% disp(['BRW: ',num2str(brw_index)]);
% disp(['WV: ',num2str(wv_index)]);
% disp(['CI: ',num2str(ci_)]);

end