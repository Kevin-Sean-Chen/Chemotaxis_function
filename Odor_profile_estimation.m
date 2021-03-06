%Odor_profile_estimation
%%% Numerical estimation of odour gradient

%%
%Specifying parameters
nx = 15;                          %Number of steps in space(x)
ny = 15;                          %Number of steps in space(y)->in a 90mm plate
nz = 15;                          %Number of steps in space(z)->10mm depth of plate
nt = 200;                         %Number of time steps 
dt = 0.01;                        %Width of each time step
dx = 1/(nx-1);                    %Width of space step(x)
dy = 1/(ny-1);                    %Width of space step(y)
dz = 1/(nz-1);                    %Width of space step(z)
x = (0:dx:1)*9;                   %Range of x(0,9) and specifying the grid points
y = (0:dy:1)*9;                   %Range of y(0,9) and specifying the grid points
z = (0:dz:1)*1;                   %Range of z(0,1) and specifying the grid points
u = zeros(nx,ny,nz);              %Preallocating u
un = zeros(nx,ny,nz);             %Preallocating un
vis = 0.081*1;                    %Diffusion coefficient/viscocity   %%%MEK diffiusion coefficient in air
% UW=0;                             %x=0 Dirichlet B.C 
% UE=0;                             %x=L Dirichlet B.C 
% US=0;                             %y=0 Dirichlet B.C 
% UN=0;                             %y=L Dirichlet B.C 
UnW=0.;                            %x=0 Neumann B.C (du/dn=UnW) west
UnE=0;                            %x=L Neumann B.C (du/dn=UnE) east
UnS=0;                            %y=0 Neumann B.C (du/dn=UnS) south
UnN=0;                            %y=L Neumann B.C (du/dn=UnN) north
UnU=0;                            %z=0 Neumann B.C (du/dn=UnU) up
UnD=0;                            %z=L Neumann B.C (du/dn=UnD) down
%%
%Initial Conditions
C0 = 2;  %initial concentration
X0 = [4.,5.,0.0,0.25];  %butanone plug
% for i=1:nx
%     for j=1:ny
%         for k=1:nz
%             if ((X0(1)<=y(j)) & (y(j)<=X0(2)) & ...
%                 (X0(1)<=x(i)) & (x(i)<=X0(2)) & ...
%                 (X0(3)<=z(k)) & (z(k)<=X0(4)) )
%                 u(i,j,k)=C0;
%             else
%                 u(i,j,k)=0;
%             end
%         end
%     end
% end
xx = find(X0(1)<=x & X0(2)>=x);
yy = find(X0(1)<=y & X0(2)>=y);
zz = find(X0(3)<=z & X0(4)>=z);
u(xx,yy,zz) = C0;

%%
%B.C vector
bc=zeros(nx-2,ny-2,nz-2);
% bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
% bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs
bc(1,:,:)=-UnW/dx; bc(nx-2,:,:)=UnE/dx;  %Neumann B.Cs
bc(:,1,:)=-UnS/dy; bc(:,ny-2,:)=UnN/dy;  %Neumann B.Cs
bc(:,:,1)=-UnU/dz; bc(:,:,nz-2)=UnD/dz;  %Neumann B.Cs
%B.Cs at the corners:
% bc(1,1)=UW/dx^2+US/dy^2;    bc(nx-2,1)=UE/dx^2+US/dy^2;
% bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
% bc(1,1,1)=0;  bc(nx-2,1,1)=0;  bc(1,nx-2,1)=0; bc(1,1,nx-2)=0;
% bc(nx-2,nx-2,1)=0; bc(nx-2,1,nx-2)=0; bc(1,nx-2,nx-2)=0; bc(nx-2,nx-2,nx-2)=0;
% bc=vis*dt*bc;
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2, 1:nx-3, 1, nx-2, nx-2);
%Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
Ax(1,1,1)=-1; Ax(nx-2,nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
%Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
Ay(1,1,1)=-1; Ay(ny-2,ny-2,ny-2)=-1;  %Neumann B.Cs
Ez=sparse(2:nz-2,1:nz-3,1,nz-2,nz-2);
%Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
Az(1,1,1)=-1; Az(nz-2,nz-2,nz-2)=-1;  %Neumann B.C
%A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2)+kron();
%D=speye((nx-2)*(ny-2))-vis*dt*A;

% e=ones(n,1); I=speye(n);
% L1=(alpha/(h^2))*spdiags([e -2*e e],-1:1,n,n);
% L = kron(L1,kron(I,kron(L1,I))+kron(I,kron(I,L1)));
% L = L + kron(L1,kron(I,I));
% A = kron();
%%
%Calculating the field variable for each time step for 3-D diffusion
figure;
for it=0:nt
    un=u;
    %h=surf(x,y,squeeze(u(:,:,3))','EdgeColor','none');       %plotting the field variable
    %h=surf(squeeze(u(:,:,4))','EdgeColor','none');       %plotting the field variable
    [X,Y,Z] = ndgrid(1:size(u,1), 1:size(u,2), 1:size(u,3));
    pointsize = 3;
    h = scatter3(X(:), Y(:), Z(:), pointsize, u(:));
    %h = slice(X,Y,Z,u,1,[],[],'nearest');
    shading interp
    %axis ([0 9 0 9 0 2])
    title({['3-D Diffusion with {\nu} = ',num2str(vis)];['time (\itt) = ',num2str(it*dt)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Transport property profile (u) \rightarrow')
    colormap jet
    colorbar;
    drawnow; 
    refreshdata(h)
    %Uncomment as necessary
    %Implicit method:
    %{
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    U=D\U;
    U=reshape(U,nx-2,ny-2);
    u(2:nx-1,2:ny-1)=U;
    %Boundary conditions
    %Dirichlet:
%     u(1,:)=UW;
%     u(nx,:)=UE;
%     u(:,1)=US;
%     u(:,ny)=UN;
    %Neumann:
    u(1,:)=u(2,:)-UnW*dx;
    u(nx,:)=u(nx-1,:)+UnE*dx;
    u(:,1)=u(:,2)-UnS*dy;
    u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
    %Explicit method:{
    for i=2:nx-1
        for j=2:ny-1
            for k=2:nz-1
                u(i,j,k)=un(i,j,k)+(vis*dt*(un(i+1,j,k)-2*un(i,j,k)+un(i-1,j,k))/(dx*dx))+...
                (vis*dt*(un(i,j+1,k)-2*un(i,j,k)+un(i,j-1,k))/(dy*dy))+...
            (vis*dt*(un(i,j,k+1)-2*un(i,j,k)+un(i,j,k-1))/(dz*dz));
            end
        end
    end
    %Boundary conditions
    %Dirichlet:
%     u(1,:)=UW;
%     u(nx,:)=UE;
%     u(:,1)=US;
%     u(:,ny)=UN;
    %Neumann:
%     u(1,:,:)=u(2,:,:)-UnW*dx;
%     u(nx,:,:)=u(nx-1,:,:)+UnE*dx;
%     u(:,1,:)=u(:,2,:)-UnS*dy;
%     u(:,ny,:)=u(:,ny-1,:)+UnN*dy;
%     u(:,:,1)=u(:,:,2)-UnU*dz;
%     u(:,:,nz)=u(:,:,nz-1)+UnD*dz;
    %}
    %modified boundary condition
    k = 1; %leakage
    A = -3.5*10^-11;  B = 1/3.19;  %effective evaporation
    u(1,:,:)=u(2,:,:)-k*u(2,:,:)*dx;
    u(nx,:,:)=u(nx-1,:,:)+k*u(nx-1,:,:)*dx;
    u(:,1,:)=u(:,2,:)-k*u(:,2,:)*dy;
    u(:,ny,:)=u(:,ny-1,:)+k*u(:,ny-1,:)*dy;
    u(:,:,1)=u(:,:,2)-(A+B*u(:,:,2))*dz;
    u(:,:,nz)=u(:,:,nz-1)+k*u(:,:,nz-1)*dz;
end

%%
rho=1.2;
cp=50;
%Length in thee three directions (i.e. width, depth and height)
W=100;
D=50;
H=30;
%Number of nodes in the three directions
NW=101;
ND=51;
NH=31;
%Distance between nodes in each of the three directions
w=W/(NW-1);
d=D/(ND-1);
h=H/(NH-1);
%Number of time steps and time interval
Nt=18000;
dt=1;
Total_t=Nt*dt;
%Creating matrix for the temperature
Tn=zeros(NH,NW,ND);
%Creating matrix for thermal conductivity
K=10*ones(NH,NW,ND);
K(:,1,:)=0.001;
K(:,end,:)=0.001;
K(1,:,:)=0.001;
K(end,:,:)=0.001;
K(:,:,1)=0.001;
K(:,:,end)=0.001;
%Initial conditions
Tn(:,:,:)=0;
t=0;
[X,Y,Z] = ndgrid(1:NH,1:NW,1:ND);
%March in time
for k=1:Nt
    k
    %Update time
    t=t+dt;
      %Update the temperature
      Tc=Tn;
      %Calculating new temperature in x and y and z locations
      for i=2:NW-1
          for j=2:ND-1
              for m=2:NH-1
                  Tn(m,i,j)=Tc(m,i,j)+...
                      dt*(K(m,i,j)/rho/cp)*...
                      (((Tc(m+1,i,j)-2*Tc(m,i,j)+Tc(m-1,i,j))/h/h)+...
                      ((Tc(m,i+1,j)-2*Tc(m,i,j)+Tc(m,i-1,j))/w/w)+...
                      ((Tc(m,i,j+1)-2*Tc(m,i,j)+Tc(m,i,j-1))/d/d));
              end
          end
      end
       scatter3(X(:),Y(:),Z(:),1,Tn(:),'filled')
       colorbar
       title(sprintf('Time Step:%d',k))
       drawnow
      %Applying the source term
      if (t<2*3600)
          Tn(20,80,40)= Tn(20,80,40)+(dt*100000/rho/cp);
      end
      %Boundary condition
      Tn(:,1,:)=Tn(:,2,:);
      Tn(:,end,:)=Tn(:,end-1,:);
      Tn(1,:,:)=Tn(2,:,:);
      Tn(end,:,:)=Tn(end-1,:,:);
      Tn(:,:,1)=Tn(:,:,2);
      Tn(:,:,end)=Tn(:,:,end-1);
end