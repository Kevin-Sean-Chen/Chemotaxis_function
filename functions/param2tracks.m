%%% function with output tracks and chemotaxis-index, with input as the
%%% model parameters
function [simdata, CI] = param2tracks(x, specs, v_dist)
%%%
% input parameter vector x with the mVM model parameters
% reconconstruct the kernels for the generative model
% ouput tracks along the landscape
%%%
    
    %%% loading parameters
    K1 = x(1);
    A_ = x(2);
    B_ = x(3:6);
    C_ = x(7);
    Amp = x(8);
    tau = x(9); 
    Amp_h = x(10); 
    tau_h = x(11); 
    K2 = x(12);  
    gamma = x(13);
    
    %%% load specifics for chemotaxis simulation
    M = specs.M;
    fr = specs.fr;
    cosBasis = specs.cosBasis;
    REP = specs.REP;
    T = specs.T;
    dt = specs.dt;
    
%     fr = 14/5;
%     [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(4, [0, 8], 1.3);
%     REP = 50;
%     T = floor(20*60*fr);   
% %     Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623.mat');
%     Cmap = load('/projects/LEIFER/Kevin/Data_odor_flow_equ/Landscape_low_0623_2.mat');
% %     M = Cmap.vq1;
%     M = fliplr(flipud(M)); 
%     dt = 1;
    
    %%% reconstructing kernels and parameters
    wind = length(cosBasis);
    Kwin = 1:wind;
    kappa = ((K1)^1)^1;
    kappa2 = ((K2)^1)^1;
    if kappa2<0.05; kappa2 = 0.1; end  %simplify calculation
    A = A_*1;
    Kdcp = Amp*exp(-Kwin/tau)*1;
    Kdth = Amp_h*exp(-Kwin/tau_h)*1;
    Kddc = B_ * cosBasis'*1;
    Pturn_base = C_;
    
    %%% record simulation
    clear simdata
    simdata(length(REP)) = struct();
    simdata(1).dth = []; %angle change
    simdata(1).dcp = [];  %perpendicular dC
    simdata(1).dc = [];  %tangential dC
    simdata(1).dis = [];  %displacements
    simdata(1).mask = [];  %mark for tracks
    simdata(1).xy = [];  %all track location in 2D
    simdata(1).theta = [];  %local angles
    
    %%% simlation loop
    [Mx, My] = size(M);

    for rep = 1:REP
        vm = 0.2*30/fr;  %should be adjusted with the velocity statistics~~ this is approximately 0.2mm X 
        vs = .1;
        perp_dist = 1;
        tracks = zeros(T,2);
        %%% from middle
        tracks(1,:) = [size(M,2)*1/2 + randn()*100, randn()*100 + size(M,1)*1/2]*1. + 0.*[size(M,2)*rand(), size(M,1)*5/8];%origin; %initial position
%         tracks(1,:) = [size(M,2)*1/2 + randn()*0, randn()*0 + size(M,1)*1/2]*1. + 0.*[size(M,2)*rand(), size(M,1)*5/8];

        tracks(2,:) = tracks(1,:)+randn(1,2)*vm*dt;%origin+randn(1,2)*vm*dt;
        ths = zeros(1,T);  ths(1:3) = randn(1,3)*360; %initial angle
        dxy = randn(1,2);  %change in each step

        %"recordings"
        dths = zeros(1,T); 
        dcs = zeros(1,T);
        dcps = zeros(1,T);

        %placer for filtering vectors
        dCv = zeros(1,wind);
        dCpv = zeros(1,wind);
        dthv = zeros(1,wind);
        for t = 2:T

            %%% dC in a window
            [xi, yi] = plate_bound(M, tracks(t-1,1), tracks(t-1,2));
            dCv = [ M(yi, xi) , dCv];
            dCv = dCv(1:end-1);  %remove out of window
            perp_dir = [-dxy(2) dxy(1)];
            perp_dir = perp_dir/norm(perp_dir)*perp_dist;
            %%% dCp in a window
            [xil, yil] = plate_bound(M, tracks(t-1,1)+perp_dir(1), tracks(t-1,2)+perp_dir(2));
            [xir, yir] = plate_bound(M, tracks(t-1,1)-perp_dir(1), tracks(t-1,2)-perp_dir(2));
            dCpv = [(M(yil,xil) - M(yir,xir))/ 1 ,  dCpv];
            dCpv = dCpv(1:end-1);
            %%% actions
            wv = (1*sum(Kdcp.*dCpv) + 0*1) + (vmrand(0,kappa))*180/pi;%kappa^1*randn;%length(wind)
            P_event = (A - Pturn_base) / (1+exp( -(sum(Kddc.*dCv)*1. + (sum(Kdth.*abs(dthv)*1)) *dt + 0)+0) ) + Pturn_base;%length(wind)
            P_event = P_event*5/14;  % rate per sec
            if rand < P_event*1
                beta = 1;
            else
                beta = 0;
            end

            if rand < gamma
                rt = beta*(vmrand(pi,1*kappa2)*180/pi);
            else
                rt = beta*(vmrand(0,0)*180/pi);
            end
            dth = wv+rt;
            dth = wrapToPi(dth*pi/180)*180/pi;

            %%% angular or turning history
            dthv = [dth*pi/pi , dthv]; %[dth  beta];%
            dthv = dthv(1:end-1);
            dths(t) = dth;
            dcs(t) = dCv(1);
            dcps(t) = dCpv(1);

            %%% draw velocity
            vv = vm+vs*randn;
%             vv = v_dist(randi(length(v_dist)));
%             if vv<1
%                 vv = 1;
%             end
            ths(t) = ths(t-1)+dth*dt;
            dd = [vv*sin(ths(t)*pi/180) vv*cos(ths(t)*pi/180)];
            dxy = dd';%(R)*dd';

            tracks(t,1) = tracks(t-1,1)+dxy(1)*dt;
            tracks(t,2) = tracks(t-1,2)+dxy(2)*dt;

            %%% within the odor environment
            if tracks(t,1)>size(M,2)-1 | tracks(t,1)<3 | tracks(t,2)>size(M,1)-1 | tracks(t,2)<3
                tracks_ = zeros(t,2);
                tracks_(:,1) = tracks(1:t,1);
                tracks_(:,2) = tracks(1:t,2);  %ending track
                tracks = tracks_;
                break;
            end
        end
        
        pos = find(dths==0 | dcs==0 | dcps==0);
        dths(pos) = [];
        dcps(pos) = [];
        dcs(pos) = [];
        
        % store as structure
        simdata(rep).dth = dths; %angle change
        simdata(rep).dcp = dcps;  %perpendicular dC
        simdata(rep).dc = dcs;  %tangential dC
        simdata(rep).xy = tracks';  %all track location in 2D
        simdata(rep).theta = ths;  %local angles
    end
    
    %%% CI
    n_up = 0;
    for rep = 1:REP
        dc = simdata(rep).dc;
        if dc(1) < dc(end)
%         if dc(end)>max(max(M))*0.5
            n_up = n_up + 1;
        end
    end
    CI = 2*n_up / REP - 1;
end