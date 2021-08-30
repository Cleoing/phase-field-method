$ git config --global core.autocrlf true
# Configure Git to ensure line endings in files you checkout are correct for Windows.
# For compatibility, line endings are converted to Unix style when you commit files.

clear;clc;clf;tic

%-----------basic set up------------
% set up area and grid size
nx=9;    
ny=nx;
dx=0.03;
xnode=dx:dx:nx;
r=length(xnode);
dy=dx;
ynode=dy:dy:ny;
s=length(ynode);

% set up time step
tf=0.2;
dt=0.0001;
tvec=dt:dt:tf; 
nt=length(tvec);

% parameters
tau=0.00025; % relaxation time
K=2.0; % lantent heat
delta=0.02; % Anisotropic strength
jq=6; % Anisotropic modulus
alpha=0.9;% angle
gamma=10; % materials parameter
chi=dt/tau; % just for simplification
epsilonbar=0.011; % interface width average
R=2; % initial radius

%----------zero matrix setup---------
P=zeros(r,s);
T=zeros(r,s);
Pn=zeros(r,s);
Tn=zeros(r,s);
Te=zeros(r,s);
m=zeros(r,s);
theta=zeros(r,s);
theta_x=zeros(r,s);
theta_y=zeros(r,s);
theta0=zeros(r,s);
epsilon=zeros(r,s);
% derivative 
P_x=zeros(r,s);
P_y=zeros(r,s);
P_xx=zeros(r,s);
P_yy=zeros(r,s);
P_xy=zeros(r,s);
% anisotropic
epsilonsqu_x=zeros(r,s);
epsilonsqu_y=zeros(r,s);
en=zeros(r,s);
eta2dthe=zeros(r,s);
etathe=zeros(r,s);
etaprimethe=zeros(r,s);
eta=zeros(r,s);
etaprime=zeros(r,s);

% ------------initial condition--------------

for i=1:r
    for j=1:s
        P(i,j)=0;
        T(i,j)=0;  % T=Treal/Tm

        Te(i,j)=1.2;
        
        if (i-r/2)^2+(j-s/2)^2<R^2
            P(i,j)=1;
            T(i,j)=1;

            
        end
        
    end
end

%------------phase field discrete----------

for t=1:nt
    t
    for i=1:r
        for j=1:s
            if i==1 && j==1
                P_x(i,j)=0;
                P_y(i,j)=0;
                P_xx(i,j)=0;
                P_yy(i,j)=0;
                
            elseif i==r && j==1
                P_x(i,j)=0;
                P_y(i,j)=0;
                P_xx(i,j)=0;
                P_yy(i,j)=0;

            elseif i==1 && j==s
                P_x(i,j)=0;
                P_y(i,j)=0;
                P_xx(i,j)=0;
                P_yy(i,j)=0;

            elseif i==r && j==s
                P_x(i,j)=0;
                P_y(i,j)=0;
                P_xx(i,j)=0;
                P_yy(i,j)=0;

            elseif i<r && i>1 && j==1
                P_x(i,j)=0.5*(P(i+1,j)-P(i-1,j))/dx;
                P_y(i,j)=0;
                P_xx(i,j)=(P(i+1,j)+P(i-1,j)-2*P(i,j))/(dx)^2;
                P_yy(i,j)=0;
 
            elseif i<r && i>1 && j==s
                P_x(i,j)=0.5*(P(i+1,j)-P(i-1,j))/dx;
                P_y(i,j)=0;
                P_xx(i,j)=(P(i+1,j)+P(i-1,j)-2*P(i,j))/(dx)^2;
                P_yy(i,j)=0;
 
            elseif j<r && j>1 && i==1
                P_x(i,j)=0;
                P_y(i,j)=0.5*(P(i,j+1)-P(i,j-1))/dy;
                P_xx(i,j)=0;
                P_yy(i,j)=(P(i,j+1)+P(i,j-1)-2*P(i,j))/(dy)^2;

            elseif j<r && j>1 && i==r
                P_x(i,j)=0;
                P_y(i,j)=0.5*(P(i,j+1)-P(i,j-1))/dy;
                P_xx(i,j)=0;
                P_yy(i,j)=(P(i,j+1)+P(i,j-1)-2*P(i,j))/(dy)^2;
  
            else
                P_x(i,j)=0.5*(P(i+1,j)-P(i-1,j))/dx;
                P_y(i,j)=0.5*(P(i,j+1)-P(i,j-1))/dy;
                P_xx(i,j)=(P(i+1,j)+P(i-1,j)-2*P(i,j))/(dx)^2;
                P_yy(i,j)=(P(i,j+1)+P(i,j-1)-2*P(i,j))/(dy)^2;
                P_xy(i,j)=(P(i+1,j+1)-P(i-1,j+1)-P(i+1,j-1)+P(i-1,j-1))/(4*dx^2);
            end
        end
    end
    % --------------Anisotropic-----------------
    for i=1:r
        for j=1:s
            theta(i,j)=atan(P_y(i,j)/(P_x(i,j)+1e-80)); % add one small constance for preventing denominator equals to zero
            eta2dthe(i,j)=-delta*jq*(2*sin(jq*theta(i,j))+delta*sin(2*jq*theta(i,j))); % d eta^2/d theta
            
            etathe(i,j)=-delta*jq*sin(jq*theta(i,j));                  % d eta/d theta
            etaprimethe(i,j)=-delta*jq^2*cos(jq*theta(i,j));            % d etaprime/d theta
            
            eta(i,j)=1+delta*cos(jq*theta(i,j));                   % eta
            
           theta_x(i,j)=(P_x(i,j)*P_xy(i,j)-P_y(i,j)*P_xx(i,j))/(P_x(i,j)^2+P_y(i,j)^2+1e-40); % d theta/dx
           theta_y(i,j)=(P_x(i,j)*P_yy(i,j)-P_y(i,j)*P_xy(i,j))/(P_x(i,j)^2+P_y(i,j)^2+1e-40); % d theta/dy
           %epsilon_d(i,j)=-epsilonbar*delta*j*sin(jq*(theta(i,j)-theta0(i,j)));
        end
    end

    epsilon=epsilonbar*eta;
    
    %-------------simulation-------------
    for i=1:r
        for j=1:s
           m(i,j)=alpha/pi*atan(gamma*(Te(i,j)-T(i,j))); 
           
           A1=P(i,j)*(1-P(i,j))*(P(i,j)-0.5+m(i,j)); % double well energy
           A2=epsilonbar^2*eta2dthe(i,j)*(P_x(i,j)*theta_x(i,j)+P_y(i,j)*theta_y(i,j));% partial for anisotropy
           A3=epsilon(i,j)^2*(P_xx(i,j)+P_yy(i,j));% partial for phase field
           A4=epsilonbar^2*(P_x(i,j)*theta_y(i,j)-P_y(i,j)*theta_x(i,j))*(etathe(i,j)^2+eta(i,j)*etaprimethe(i,j));% angle
           
           Pn(i,j)=P(i,j)+chi*(A1+A2+A3+A4);
           
  % --------------time discrete--------------         
           if i==1 && j==1
                Tn(i,j)=K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif i==r && j==1
                Tn(i,j)=K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif i==1 && j==s
                Tn(i,j)=K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif i==r && j==s
                Tn(i,j)=K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif i<r && i>1 && j==1
                Tn(i,j)=dt/(dx)^2*(T(i+1,j)+T(i-1,j)-2*T(i,j))+K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif i<r && i>1 && j==s
                Tn(i,j)=dt/(dx)^2*(T(i+1,j)+T(i-1,j)-2*T(i,j))+K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif j<r && j>1 && i==1
                Tn(i,j)=dt/(dy)^2*(T(i,j+1)+T(i,j-1)-2*T(i,j))+K*(Pn(i,j)-P(i,j))+T(i,j);
            elseif j<r && j>1 && i==r
                Tn(i,j)=dt/(dy)^2*(T(i,j+1)+T(i,j-1)-2*T(i,j))+K*(Pn(i,j)-P(i,j))+T(i,j);
           else
                 Tn(i,j)=dt*(T(i,j+1)+T(i,j-1)+T(i+1,j)+T(i-1,j)+0.5*(T(i+1,j+1)+T(i+1,j-1)+T(i-1,j+1)+T(i-1,j-1))-6*T(i,j))/(dx)^2+K*(Pn(i,j)-P(i,j))+T(i,j);
            end
            
        end
    end
    P=Pn;
    T=Tn; 
end
%------------figure--------------
x=1:r;
y=1:s;
t=nt;
Ppost=P;




contourf(x,y,Ppost)

axis square
