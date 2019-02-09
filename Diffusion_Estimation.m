%Diffusion estimation

% Create a model.
model = createpde(1);
% Create geometry and assign it to the model
C1 = [1;-1;0;.6];
%C2 = [1;-1;0;0.5];  %circle==1; x-y cneter; radius
%C2 = [C2;zeros(length(R1) - length(C2),1)];
gd = [C1,C1];%[C1,C2];
sf = 'R1+C1';
ns = char('R1','C1');%,'C1')';
g = decsg(gd,sf,ns);
geometryFromEdges(model,g);
figure
pdegplot(model,'EdgeLabels','on','FaceLabels','on')
% Specify coefficients, diffusion coefficient as c and source as f
% Assign zero source and a diffusion coefficient, c = 1
specifyCoefficients(model,'Face',1,'m',0,'d',1,'c',1,'a',0,'f',0);
% Assign non-zero source and a different diffusion coefficient c = 2 on inner
% circle 
specifyCoefficients(model,'Face',2,'m',0,'d',1,'c',2,'a',0,'f',1);
% Set initial condition
setInitialConditions(model,0);
setInitialConditions(model,1,'Face',2);
% Generate mesh and solve.
generateMesh(model)
tlist = linspace(0,0.1,20);
result = solvepde(model,tlist)
% Plot results of last time-step.
figure
pdeplot(model,'XYData',result.NodalSolution(:,end))

%% 
% Simulating the 2-D Diffusion equation by the Finite Difference
...Method 
% Numerical scheme used is a first order upwind in time and a second
...order central difference in space (Implicit and Explicit)
%%
%Specifying parameters
nx=100;                          %Number of steps in space(x)
ny=100;                          %Number of steps in space(y)       
nt=100;                           %Number of time steps 
dt=0.01;                         %Width of each time step
dx=9/(nx-1);                     %Width of space step(x)  ###9cm plate
dy=9/(ny-1);                     %Width of space step(y)
x=0:dx:9;                        %Range of x(0,10) and specifying the grid points
y=0:dy:9;                        %Range of y(0,10) and specifying the grid points
u=zeros(nx,ny);                  %Preallocating u
un=zeros(nx,ny);                 %Preallocating un
vis=0.081;                       %Diffusion coefficient/viscocity   %%%MEK diffiusion coefficient in air
UW=0;                            %x=0 Dirichlet B.C 
UE=0;                            %x=L Dirichlet B.C 
US=0;                            %y=0 Dirichlet B.C 
UN=0;                            %y=L Dirichlet B.C 
UnW=0;                           %x=0 Neumann B.C (du/dn=UnW)
UnE=0;                           %x=L Neumann B.C (du/dn=UnE)
UnS=0;                           %y=0 Neumann B.C (du/dn=UnS)
UnN=0;                           %y=L Neumann B.C (du/dn=UnN)
%%
%Initial Conditions
C0 = 2;  %initial concentration
X0 = [4.3,4.7];  %butanone plug
for i=1:nx
    for j=1:ny
        if ((X0(1)<=y(j))&&(y(j)<=X0(2))&&(X0(1)<=x(i))&&(x(i)<=X0(2)))
            u(i,j)=C0;
        else
            u(i,j)=0;
        end
    end
end
%%
%B.C vector
bc=zeros(nx-2,ny-2);
bc(1,:)=UW/dx^2; bc(nx-2,:)=UE/dx^2;  %Dirichlet B.Cs
bc(:,1)=US/dy^2; bc(:,ny-2)=UN/dy^2;  %Dirichlet B.Cs
%bc(1,:)=-UnW/dx; bc(nx-2,:)=UnE/dx;  %Neumann B.Cs
%bc(:,1)=-UnS/dy; bc(:,nx-2)=UnN/dy;  %Neumann B.Cs
%B.Cs at the corners:
bc(1,1)=UW/dx^2+US/dy^2; bc(nx-2,1)=UE/dx^2+US/dy^2;
bc(1,ny-2)=UW/dx^2+UN/dy^2; bc(nx-2,ny-2)=UE/dx^2+UN/dy^2;
bc=vis*dt*bc;
%Calculating the coefficient matrix for the implicit scheme
Ex=sparse(2:nx-2,1:nx-3,1,nx-2,nx-2);
Ax=Ex+Ex'-2*speye(nx-2);        %Dirichlet B.Cs
%Ax(1,1)=-1; Ax(nx-2,nx-2)=-1;  %Neumann B.Cs
Ey=sparse(2:ny-2,1:ny-3,1,ny-2,ny-2);
Ay=Ey+Ey'-2*speye(ny-2);        %Dirichlet B.Cs
%Ay(1,1)=-1; Ay(ny-2,ny-2)=-1;  %Neumann B.Cs
A=kron(Ay/dy^2,speye(nx-2))+kron(speye(ny-2),Ax/dx^2);
D=speye((nx-2)*(ny-2))-vis*dt*A;
%%
%Calculating the field variable for each time step
i=2:nx-1;
j=2:ny-1;
for it=0:nt
    un=u;
    h=surf(x,y,u','EdgeColor','none');       %plotting the field variable
    shading interp
    axis ([0 9 0 9 0 2])
    title({['2-D Diffusion with {\nu} = ',num2str(vis)];['time (\itt) = ',num2str(it*dt)]})
    xlabel('Spatial co-ordinate (x) \rightarrow')
    ylabel('{\leftarrow} Spatial co-ordinate (y)')
    zlabel('Transport property profile (u) \rightarrow')
    drawnow; 
    refreshdata(h)
    %Uncomment as necessary
    %Implicit method:
    U=un;U(1,:)=[];U(end,:)=[];U(:,1)=[];U(:,end)=[];
    U=reshape(U+bc,[],1);
    U=D\U;
    U=reshape(U,nx-2,ny-2);
    u(2:nx-1,2:ny-1)=U;
    %Boundary conditions
    %Dirichlet:
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    %Neumann:
    %u(1,:)=u(2,:)-UnW*dx;
    %u(nx,:)=u(nx-1,:)+UnE*dx;
    %u(:,1)=u(:,2)-UnS*dy;
    %u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
    %Explicit method:
    %{
    u(i,j)=un(i,j)+(vis*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+(vis*dt*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(dy*dy));
    %Boundary conditions
    %Dirichlet:
    u(1,:)=UW;
    u(nx,:)=UE;
    u(:,1)=US;
    u(:,ny)=UN;
    %Neumann:
    %u(1,:)=u(2,:)-UnW*dx;
    %u(nx,:)=u(nx-1,:)+UnE*dx;
    %u(:,1)=u(:,2)-UnS*dy;
    %u(:,ny)=u(:,ny-1)+UnN*dy;
    %}
end