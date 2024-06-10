% BR_analysis
%%%
% analyze effective velocity and Brownian ratchet
% sample tracks and compare strains
%%%
rng(37) %37, 13

%% load data
%%% for strains
% load('/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/Data_nai_staPAW.mat') %%% odor naive
% load('/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/Data_ave_staPAW.mat') %%% odor aversive
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIZ_ave_staPAW_vars.mat') %4500
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIB_app_staPAW_vars.mat')  %200000
% load('/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIY_app_staPAW_vars.mat') %130000

strain_data = {'/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/Data_nai_staPAW.mat',...
               '/projects/LEIFER/Kevin/Data_learn/N2/ssm_analysis/Data_ave_staPAW.mat',...
               '/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIB_nai_staPAW_vars.mat',...
               '/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIY_nai_staPAW_vars.mat',...
               '/projects/LEIFER/Kevin/Publications/Chen_states_2024/Data_AIZ_ave_staPAW_vars.mat'};  % naive, ave, AIB, AIY, AIZ
%%% for flow
Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
M = Cmap.vq1;
M = fliplr(flipud(M));

%% parameter settings
pre_t = 20; %8, 2, 12, 20
post_t = 20;
n_repeat = 30;  % repeat samples
n_samp = 30;  % number of tracks
n_strains = length(strain_data);
target = [1500  2500];   % odor source
scale2speed = (1/30)/(5/14);

%% looping analysis (across strains, repeats; for CI and BR analysis)
%%% measurements
v_ci = zeros(n_strains,n_repeat);
v_br = zeros(n_strains,n_repeat);

for ww = 1:n_strains
    
    ww  %% across strains
    load(strain_data{ww});   % load strains!
    n_tracks = length(Data);

for ss = 1:n_repeat   %%%%%%% resample loop
    
%         ss %%% repeats
        temp = randperm(n_tracks);
        Data_i = Data(temp(1:n_samp));  % sub-sampling tracks
        [xxi, yyi, triali, timei] = data2xy(Data_i);
        disi = extractfield(Data_i, 'dis');  % Data
        yyi = [yyi; disi];
        maski = triali;
        [logp_test,gams_,xis_test,xisum_test,logcs] = runFB_GLMHMM_xi(mmhat,xxi,yyi,maski);
        xyi = [];
        for ii = 1:length(Data_i);  xyi = [xyi  Data_i(ii).xy]; end
        
        %%% asigning states
        S_state = argmax(sum(gams_,2));
        T_state = argmin(sum(gams_,2));
        
        %%% compute bearing
        state_vec = gams_(T_state,:)*0;
%         if ww==3; state_vec = gams_(2,:)*0;; end  % dealing with flipping
        pos_state_stapaw = find(gams_(1,:)>0.5); %<gams_(1,:));  %%% check this
        pos_state = pos_state_stapaw; 
        state_vec(pos_state) = ones(1,length(pos_state));
        trans_pos = diff(state_vec);
        trans12 = find(trans_pos>0)-0;
        trans21 = find(trans_pos<0)+0;
        npairs = min([length(trans12), length(trans21)]);
        Bpairs = zeros(2,npairs);
        vecs = diff(xyi,1,2);
        for bb = 8:npairs-8
            %%% pre vector
            v1 = (target - xyi(:, trans12(bb))');  %target
            v2 = -(xyi(:,trans12(bb)) - xyi(:,trans12(bb)-pre_t))';  % position-based  
            Bpre = angles(v1, v2);

            %%% post vector
            v1 = (target - xyi(:, trans21(bb))'); %target
            v2 = -(xyi(:,trans21(bb)) - xyi(:,trans21(bb)+post_t))';
            Bpost = angles(v1, v2);

            %%% recording pairs
            Bpairs(:,bb) = [Bpre; Bpost];
        end
        
        %%% compute alignment probability
        test = Bpairs(2,:);
        P_plus = length(find(abs(test)<90))/length(test);
        test = Bpairs(1,:);
        P_minus = length(find(abs(test)<90))/length(test);

        %%% compute run speed
        pos_run = find(gams_(S_state,:)>0.5); %<gams_(1,:));  %%% check this
        LT = mean(yyi(2,pos_run));
        frac_run = length(pos_run)/length(yyi);
        
        BR_i = LT*(P_plus - P_minus)*(frac_run)* (scale2speed);  % effective drift in Brownian ratchet
        v_br(ww, ss) = BR_i*1;
        
        %%% comput CI for tracks
        n_up = 0;  n_down = 0;
        for ii = 1:n_samp
            xyii = Data_i(ii).xy;
            dC = -(M(floor(xyii(2,1)),floor(xyii(1,1))) - M(floor(xyii(2,end)),floor(xyii(1,end))));  % by C 
%             dC =  -( (sum(([floor(xyii(2,1)),floor(xyii(1,1))] - [target]).^2))^0.5 - (sum(([floor(xyii(2,end)),floor(xyii(1,end))] - [target]).^2))^0.5 ); % by location
            if dC>=0
                n_up = n_up+1;%
%                 n_up = n_up + (sum(([floor(xyii(2,1)),floor(xyii(1,1))] - [floor(xyii(2,end)),floor(xyii(1,end))]).^2))^0.5;
            else
                n_down = n_down+1;
%                 n_down = n_down + (sum(([floor(xyii(2,1)),floor(xyii(1,1))] - [floor(xyii(2,end)),floor(xyii(1,end))]).^2))^0.5;
            end
        end
        ci_i = (n_up-n_down) / (n_up+n_down);  % sample chemotaxis index
        v_ci(ww, ss) = ci_i*1;
                
end

end
%% plotting
figure()
labels = {'N2 naive', 'N2 aversive', 'AIB-', 'AIY-', 'AIZ-'};
cols = ['k','r','b','g','c'];
for ww = 1:n_strains
    plot(v_br(ww,:), v_ci(ww,:), '.','MarkerSize',15, 'color',cols(ww)); hold on
end
legend(labels);
xlabel('v_{drift} (Brownian Ratchet)'); ylabel('CI')
set(gcf,'color','w'); set(gca,'Fontsize',20);

%%% regression
X = [v_br(:), ones(length(v_br(:)),1)];
Y = v_ci(:);
beta = X\Y;
yCalc1 = X*beta;
hold on
[~,bb]=sort(v_br(:));
plot(sort(v_br(bb)),yCalc1(bb), 'k--','LineWidth',2)
Rsq1 = 1 - sum((Y - yCalc1).^2)/sum((Y - mean(Y)).^2)
