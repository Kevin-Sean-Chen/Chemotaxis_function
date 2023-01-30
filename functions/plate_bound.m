function [xi, yi] = plate_bound(M, x, y)
%%%
% Given the odor map space M matrix and the track position x
% return the integar index in matrix M that is positive and bounded by M
%%%

[My,Mx] = size(M);  %important that the matrix has a transposed x-y!
xi = min(max(floor(x),1),Mx);
yi = min(max(floor(y),1),My);

end