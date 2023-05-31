function [xx,yy,mask, time] = data2xy(Data)
%%%
% Input Data structure that contains all tracks and the corresponding
% fields. This function extracts the relevant fields and concatenates it
% into a long vector for group-level model fitting.

% Input: Data structure
% Output: xx as the input vector [2 x N*T], with N tracks and T time
%         yy as the output vector [1 x N*T]
%         mask as a logic vector to mask transitions
%%%

% extract fields
allas = extractfield(Data, 'dth');
alldC = extractfield(Data, 'dc');
alldcp = extractfield(Data, 'dcp');
alltrials = extractfield(Data, 'mask');
time = extractfield(Data, 'time');

% make vectors
yy = allas;
xx = [alldC; 
      alldcp];
mask = true(1,length(yy));
mask(isnan(alltrials)) = false;

end