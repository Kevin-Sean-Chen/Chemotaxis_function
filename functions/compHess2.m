function [H, g] = compHess2(fun, x0, dx, Data)
%  [H, g] = compHess(fun, x0, dx, varargin)
%
%  Numerically computes the Hessian of a function fun around point x0
%  expects fun to have sytax:  y = fun(x, varargin);

n = length(x0);
H = zeros(n,n);
g = zeros(n,1);
f0 = fun(x0, Data);
A = diag(dx*ones(n,1)/2);

for j = 1:n  % compute diagonal terms
    f1 = fun(x0+2*A(:,j),Data);
    f2 = fun(x0-2*A(:,j),Data);
    H(j,j) = f1+f2-2*f0;
    g(j) = (f1-f2)/2;
end

fprintf('Computing Hessian:\nrow ');
for j = 1:n-1       % compute cross terms
    fprintf(1, ' %d',j);
    for i = j+1:n
        f11 = fun(x0+A(:,j)+A(:,i),Data);
        f22 = fun(x0-A(:,j)-A(:,i),Data);
        f12 = fun(x0+A(:,j)-A(:,i),Data);
        f21 = fun(x0-A(:,j)+A(:,i),Data);
        
        H(j,i) = f11+f22-f12-f21;
        H(i,j) = H(j,i);
    end
end
fprintf('\n');

H = H/(dx.^2);
g = g/dx;