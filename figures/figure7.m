% Figure 7
% USER INSTRUCTION: download and unzip the demo data pack Chen_learning_2023.zip from figshare to a local directory of your choice, and modify datadir below accordingly. If you don't modify datadir, this script assumes it's in your system's default Downloads directory.
% datadir = fullfile(getenv('USERPROFILE'),'Downloads','Chen_learn_2023');
datadir = fullfile('/projects/LEIFER/Kevin/Publications/','Chen_learning_2023');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 7
% extracted from script 'mutant_fit.m' and 'mutant_var.m'

%% load full data
load(fullfile(datadir,'data4plots', 'mutant_fit_vars_std.mat'));  % containing all required fields of data

%%
%% arrow plots
figure
N2_param = squeeze(mle_params(2,:,:));  % loading N2 fitted parameters
N2_param = N2_param([1,3,2],:);
arrowHeadSize = 2;
ii = 3*5;  % all conditions
jj = 1;
test = zeros(3,5);
for cc = 1:3
    for ss = 1:5
        %%% load data
        load(data_mut{cc,ss});
        Data_n2_i = load(data_n2{cc});
        
        %%% plot arrow
        temp_param = squeeze(mut_param(cc,ss,:))';
        temp_n2_param = N2_param(cc,:);
        temp_n2_ci = N2_ci(cc);
        
        %%% convolution
        d_brw = BRW_Kc_conv(temp_param, Data) - BRW_Kc_conv(temp_n2_param, Data_n2_i.Data); 
        d_wv =  WV_Kcp_conv(temp_param, Data) - WV_Kcp_conv(temp_n2_param, Data_n2_i.Data); 
        d_ci = CIs(cc,ss) - temp_n2_ci;
        test(cc,ss) = d_brw;
        
        %%% N2 measurements
        Data_n2_i = load(data_n2{cc});
        brw_n2 = BRW_Kc_conv(temp_n2_param, Data_n2_i.Data);
        wv_n2 = WV_Kcp_conv(temp_n2_param, Data_n2_i.Data);
        
        subplot(131)
        plot(N2_ci(cc), ii,'o'); hold on
        h0 = quiver(temp_n2_ci, ii, d_ci, 0,0,  'autoScale', 'off'); hold on
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('CI'); yticks([])
        subplot(132)
        semilogx(brw_n2, ii,'o'); hold on
%         plot(brw_n2, ii,'o'); hold on
        h1 = quiver(brw_n2, ii, (d_brw), 0,0,  'autoScale', 'off'); hold on
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('BRW'); yticks([])
        subplot(133)
        plot(wv_n2, ii,'o'); hold on
        h2 = quiver(wv_n2, ii, d_wv, 0,0,  'autoScale', 'off'); hold on
        set(gcf,'color','w'); set(gca,'Fontsize',20); title('WV'); yticks([])
        ii = ii-1;
        
    end    
end

%% linear contribution analysis
P = ones(n_strains, n_strains);  % neuron x neuron perturbation matrix
P(1:size(P,1)+1:end) = 0;
B = zeros(2, n_strains);  % behavior x neuron
Ws = zeros(3,2, n_strains);  % condition x behavior x neuron
xtickLabels = {'AIA','AIB','AIY','AIZ','RIA'};
ytickLabels = {'BRW','WV'};
ttls = {'appetitive', 'aversive', 'naive'};

%%% color scheme
numColors = 256;  % Number of colors in the colormap
r = [linspace(1, 1, numColors), linspace(1, 0.5, numColors)];  % Red component
g = [linspace(0, 1, numColors), linspace(0.5, 1, numColors)];  % Green component
b = [linspace(0, 1, numColors), linspace(1, 1, numColors)];  % Blue component
customMap = [r', g', b'];

figure
for cc = 1:3
    for ss = 1:5    
        %%% normalized index
        load(data_mut{cc,ss});
        temp_param = squeeze(mut_param(cc,ss,:))';
        B(1,ss) = BRW_Kc_conv(temp_param, Data);
        B(2,ss) = WV_Kcp_conv(temp_param, Data);
    end
    B = B./vecnorm(B,1,2);  % normalize behavior readout?
    wtemp = B/P;%B*inv(P);
    Ws(cc,:,:) = wtemp;
    
    subplot(1,3,cc)
    imagesc(wtemp); xticklabels(xtickLabels); yticklabels(ytickLabels); xticks([1:n_strains]); yticks([1:2]); title(ttls(cc)); caxis([-0.3,.3])
    set(gcf,'color','w'); set(gca,'Fontsize',20);
%     caxis([-max(abs(wtemp(:))), max(abs(wtemp(:)))]);
end

%% functions required
function [varargout] = BRW_Kc_conv(x, Data);
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    tt = 0:length(cosBasis)-1;
    Kc = x(3:6)*cosBasis';
    [xx_train, yy_train, mask_train] = data2xy(Data);%(Data(1:400));%(Data(idvec(1:400)));
    ddc_fit = xx_train(1,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
%     brw = nanstd(conv_kernel(ddc_fit.*trials_fit, Kc));
    K_h_rec = x(10)*exp(-tt/x(11));
    filt_dth = conv_kernel(abs(yy_train), K_h_rec);
    filt_dc = conv_kernel(ddc_fit.*trials_fit, Kc/1);
    dc_dth = filt_dc+filt_dth*1;
    compute_nl = (x(2)-x(7))./(1+exp(-dc_dth - 0*(filt_dth))) + x(7);
    test_nl = (x(2)-x(7))./(1+exp(-ddc_fit - 0*(filt_dth))) + x(7);
    if nargout > 1
        varargout{1} = norm(Kc)/nanstd(ddc_fit.*1);%norm(Kc);   %(nanstd(compute_nl));% - nanstd(filt_dc);
        temp = (conv_kernel(ddc_fit.*trials_fit, Kc));
        varargout{2} = dc_dth;
    else
        varargout{1} = norm(Kc)/nanstd(ddc_fit.*1);%norm(Kc);   %(nanstd(compute_nl));% - nanstd(filt_dc);
    end
end

function [varargout] = WV_Kcp_conv(x, Data)
    [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
    tt = 0:length(cosBasis)-1;
    Kcp = x(8)*exp(-tt/x(9));
    [xx_train, yy_train, mask_train] = data2xy(Data);%(Data(1:400));%(Data(idvec(1:400)));
    dcp_fit = xx_train(2,:);
    trials_fit = ones(1,length(mask_train));
    trials_fit(find(mask_train==0)) = NaN;
%     wv = nanstd(conv_kernel(dcp_fit.*trials_fit, Kcp));
    if nargout > 1
        varargout{1} = nanstd(conv_kernel(dcp_fit.*trials_fit, Kcp));
        temp = (conv_kernel(dcp_fit.*trials_fit, Kcp));
        varargout{2} = temp(~isnan(temp));
    else
        varargout{1} = nanstd(conv_kernel(dcp_fit.*trials_fit, Kcp));
    end
end
