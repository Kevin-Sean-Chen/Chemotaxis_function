function [v_orthogonal, v_parallel, speed, curve] = track2speed(x,y)
%%%
% input x and y coordinates, both with 1 x T vector
% output velocity along the orthogonal and parallel vectors, speed and
% curvature calculated from this decompistion, all with 1 x (T-2) vector
%%%

dx = diff(x);
dy = diff(y);

v_orthogonal = (dy(1:end-1) .* dx(2:end) - dy(2:end).* dx(1:end-1)) ./ sqrt( dx(1:end-1).^2 + dy(1:end-1).^2 );
v_parallel = (dx(1:end-1) .* dx(2:end) - dy(2:end).* dy(1:end-1)) ./ sqrt( dx(1:end-1).^2 + dy(1:end-1).^2 );
speed = sqrt( v_orthogonal.^2 + v_parallel.^2 );
curve = v_orthogonal ./ ( dx(1:end-1).^2 + dy(1:end-1).^2 );


end