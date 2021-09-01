% test for phase field discrete
% to set R=10, larger than setup radius
% let phase field parameter range from 0 to 1
% test if the discrete euqation makes sense
clear;clc;clf;tic
nx=9;    
ny=nx;
dx=0.03;
xnode=dx:dx:nx;
r=length(xnode);
dy=dx;
ynode=dy:dy:ny;
s=length(ynode);

R=10;

P=ones(r,s);
P_x=zeros(r,s);
P_y=zeros(r,s);
P_xx=zeros(r,s);
P_yy=zeros(r,s);
P_xy=zeros(r,s);


for n=0.1:0.1:1
    n
    for i=1:r
    for j=1:s
        P(i,j)=0;
        if (i-r/2)^2+(j-s/2)^2<R^2
            P(i,j)=n;    
        end
    end
    end
end






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