%Convection_Diffusion

D = 1;  %diffusion coefficient (cm^2/s)
V = 2;  %~~cm/s (?)
C0 = 100;  %mM (?)
k = 10;  %outlet

T = 500;  %seconds
dt = 0.01;
time = 0:dt:T;
nt = length(time);
dy = .5;  %cm length
yy = -10:dy:10;
ny = length(yy);
dx = .5;
xx = 0:dx:15;
nx = length(xx);
C = zeros(nx,ny,nt);
mid = round(ny/2);

C(1,[mid-2:mid+2],1) = C0;

figure;
for ti = 1:nt-1
    
    %%%Convection-Diffusion dynamics
    for xi = 2:nx-1
        for yi = 2:ny-1
            C(xi,yi,ti+1) = C(xi,yi,ti) + D*dt/(dx*dx)*(C(xi+1,yi,ti)-2*C(xi,yi,ti)+C(xi-1,yi,ti))...
                + D*dt/(dy*dy)*(C(xi,yi+1,ti)-2*C(xi,yi,ti)+C(xi,yi-1,ti))...
                - V*dt/(2*dx)*(C(xi+1,yi,ti)-C(xi-1,yi,ti))...
                - 0*dt/(2*dx)*(C(xi,yi+1,ti)-C(xi,yi-1,ti));   %no velocity along y (??)
        end
    end
    
    %%%Boundary conditions
    C(1,:,ti) = 0;  %clean air
    C(1,[mid-2:mid+2],ti+1) = C0;  %source
    Q = sum(C(1,:,ti+1));
    C(nx,:,ti) = 0;%C(nx-1,:,ti)-k*C(nx-1,:,ti)*dx;  %outlet
    C(:,1,ti) = 0;%C(:,2,ti)-k*C(:,2,ti)*dy;  %ideal boundary
    C(:,ny,ti) = 0;%C(:,ny-1,ti)-k*C(:,ny-1,ti)*dy;
    
    imagesc(squeeze(C(:,:,ti))');  pause();
    
end

%%%sptaiol section
figure;
plot(squeeze(C(:,:,ti)'))

%         k = 1; %leakage
%     A = -3.5*10^-11;  B = 1/3.19;  %effective evaporation
%     u(1,:,:)=u(2,:,:)-k*u(2,:,:)*dx;
%     u(nx,:,:)=u(nx-1,:,:)+k*u(nx-1,:,:)*dx;
%     u(:,1,:)=u(:,2,:)-k*u(:,2,:)*dy;
%     u(:,ny,:)=u(:,ny-1,:)+k*u(:,ny-1,:)*dy;
