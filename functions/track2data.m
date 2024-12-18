function [Data] = track2data(Tracks, cand, M)

Data(length(cand)) = struct();  % all tracks and chemotaxis time series

%%% pre-processing parameters
poly_degree = 3;  %polynomial fit for the tracks
bin = 5;  %temporal binning  %~0.5 s
filt = 14*1;  %filtering tracks
l_window = 1;  %lag time
perp_dist = 1;  %perpendicular vectors for C^perp measurements
fr2sec = 1/14;
dis_thr = 5.5*fr2sec*30;  %max or min mm/s that is reasonable
%basis for filters
nB = 4;
[cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, [0, 8], 1.3);

%%% gather time series
Data(1).dth = []; %angle change
Data(1).dcp = [];  %perpendicular dC
Data(1).dc = [];  %tangential dC
Data(1).dis = [];  %displacements
Data(1).mask = [];  %mark for tracks
Data(1).xy = [];  %all track location in 2D
Data(1).theta = [];  %local angles
Data(1).time = [];  % real time in the assay
Data(1).opto = [];  % optogenetics

ntracks = length(cand);  %where cand is candidates checked somewhere else
for ii = 1:ntracks
    
    %%% process path to angles and distances
    temp = Tracks(cand(ii));  %for loading saved tracks
    temp1 = zeros(round(size(temp.Path,1)/1),2);
    temp1(:,1) = smooth(temp.SmoothX, filt,'sgolay',poly_degree);
    temp1(:,2) = smooth(temp.SmoothY, filt,'sgolay',poly_degree);
    subs = temp1(1:bin:end,:);
    timev = temp.Time;
    timev = timev(1:bin:end);
    vecs = diff(subs);
%     vecs = [vecs(1,:); vecs];   % compensating for the size change/time shift after taking difference 
    angs = zeros(1,size(vecs,1));    
    ang_loc = zeros(1,size(vecs,1));
    opto_i = temp.LEDVoltages;
    opto_i = opto_i(1:bin:end);

    %%% for time step calculations
    dCp = zeros(1,size(vecs,1));
    dCs = zeros(1,size(vecs,1));
    dds = zeros(1,size(vecs,1));
    trials = ones(1,size(vecs,1));
    xys = zeros(2,size(vecs,1));
    
    %%% iterate through worms
    for dd = (l_window+1):length(angs)
        %%% angle function
        angs(dd) = angles(vecs(dd-l_window,:)/norm(vecs(dd-l_window,:)),vecs(dd,:)/norm(vecs(dd,:)));%/(norm(vecs(dd,:)));  %curving rate?
        ang_loc(dd) = angles([1,0],vecs(dd,:)/norm(vecs(dd,:)));
        
        %%% perpendicular concentration change
        perp_dir = [-vecs(dd-0,2), vecs(dd-0,1)];
        perp_dir = perp_dir/norm(perp_dir)*perp_dist;
        [xil, yil] = plate_bound(M, subs(dd-0,1)+perp_dir(1), subs(dd-0,2)+perp_dir(2));
        [xir, yir] = plate_bound(M, subs(dd-0,1)-perp_dir(1), subs(dd-0,2)-perp_dir(2));
        dCp(dd) = (M(yil,xil) - M(yir,xir));

        %%% forward concentration
        dCs(dd) = M(floor(subs(dd-0,2)), floor(subs(dd-0,1)));
        
        %%% check displacement
        dds(dd) = norm([subs(dd,1)-subs(dd-l_window,1), subs(dd,2)-subs(dd-l_window,2)]);
        
        %%% record location in 2D
        xys(:,dd) = subs(dd,:)';

    end
    %remove zeros
    dCp = dCp((l_window+1):end);
    dCs = dCs((l_window+1):end);
    dds = dds((l_window+1):end);
    angs = angs((l_window+1):end);
    ang_loc = ang_loc(l_window+1:end);
    trials = trials((l_window+1):end);
    xys = xys(:, (l_window+1):end);
    timev = timev((l_window+2):end);  % one index off because of the derivative
    opto_i = opto_i((l_window+2):end);
    
    trials(1:50) = NaN; trials(end) = NaN;
    
    % store as structure
    Data(ii).dth = angs; %angle change
    Data(ii).dcp = dCp;  %perpendicular dC
    Data(ii).dc = dCs;  %tangential dC
    Data(ii).dis = dds;  %displacements
    Data(ii).mask = trials;  %mark for tracks
    Data(ii).xy = xys;  %all track location in 2D
    Data(ii).theta = ang_loc;  %local angles
    Data(ii).Basis = cosBasis;
    Data(ii).lambda = .1;
    Data(ii).time = timev;
    Data(ii).opto = opto_i;  %opto-stimuli
    
end

end