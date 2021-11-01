%Estimate_Concentration
function [X] = Est_con(x,y,x0,y0,C0)
%location x,y
%origin x0,y0 and factor C0

%C0 = 0.2;
D = 1.5*1e-3;
duT = 60*60*100;
d = 0.18;
%2-D diffusion
%%% C = N o e?r2 /400Dt/4?dDt, 
X = C0/(4*pi*d*D*duT)*exp(-((x-x0)^2 + (y-y0)^2)/(400*D*duT));

end

%%% test
% C0 = 0.2;
% x0 = 50;
% y0 = 50;
% xx = 0:1:100;
% yy = 0:1:100;
% XX = zeros(length(xx),length(yy));
% for i = 1:length(xx)
%     for j = 1:length(yy)
%         XX(i,j) = Est_con(xx(i),yy(j),x0,y0,C0);
%     end
% end