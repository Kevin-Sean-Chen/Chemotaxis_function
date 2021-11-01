function [cosBasis, tgrid, basisPeaks] = makeRaisedCosBasis(nB, peakRange, logOffset)
% Bprs.nBasis = 8;  % number of basis vectors
% Bprs.peakRange = [0, 10]; % location of 1st and last peaks
% Bprs.dt = 0.1; % time bin size
% Bprs.logScaling = 'log';  % specify log scaling of time axis
% Bprs.logOffset = 1.5;  % nonlinear stretch factor (larger => more linear)

% Define function for single raised cosine basis function
raisedCosFun = @(x,ctr,dCtr)((cos(max(-pi,min(pi,(x-ctr)*pi/dCtr/2)))+1)/2);


% Define nonlinear time axis stretching function and its inverse
nlin = @(x)(log(x+1e-20));
invnl = @(x)(exp(x)-1e-20);

% Compute location for cosine basis centers
logPeakRange = nlin(peakRange+logOffset);   % 1t and last cosine peaks in stretched coordinates
dCtr = diff(logPeakRange)/(nB-1);   % spacing between raised cosine peaks
Bctrs = logPeakRange(1):dCtr:logPeakRange(2);  % peaks for cosine basis vectors
basisPeaks = invnl(Bctrs);  % peaks for cosine basis vectors in rescaled time

% Compute time grid points
dt = 1;
minT = 0; % minimum time bin (where first basis function starts)
maxT = invnl(logPeakRange(2)+2*dCtr)-logOffset; % maximum time bin (where last basis function stops)
tgrid = (minT:dt:maxT)'; % time grid
nT = length(tgrid);   % number of time points in basis

% Make the basis
cosBasis = raisedCosFun(repmat(nlin(tgrid+logOffset), 1, nB), repmat(Bctrs, nT, 1), dCtr);

end
