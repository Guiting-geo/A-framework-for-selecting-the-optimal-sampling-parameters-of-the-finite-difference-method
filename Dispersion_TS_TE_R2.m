clc;
clear all;

M=8;
N=8;
dx=6;
dt=0.001;

beta=0.01:0.01:pi;

theta=0.01:0.01:pi;

v_stencil=[1500,1700,2000,2500,3000];
r_stencil=[0.2,0.25,0.3,0.35, 0.4];
vt=2500;

figure;
for ir=1:5
    r=r_stencil(ir);
    rt=0.3;
    %c=fdcoeff_time_space_angles(M,2,v,dt,dx);
    c=fdcoeff_time_space_angles_r(M,0,rt)
    [m n]=size(beta);
    [nn m]=size(theta);
    delta=zeros(n,m);
    for i=1:n    %beta
        bb=beta(i);
        for j=1:m     %theta
            tt=theta(j);
            temp=c(1);
             for k=1:M
               temp=temp+2*c(1+k)*(   cos(k*bb*cos(tt)) + cos(k*bb*sin(tt))          );
             end
            delta(i,j)=1/r/bb*acos( 1+  0.5*r^2*temp              );
        end
    end
    hold on;
    plot(theta,delta(270,:),'LineWidth',1.6)     %  14/16
end

grid on;
set(gca,'linewidth',1.2);
box on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
%axis([0,pi,0.89,1.025])
axis([0,pi,0.91,1.045])
set(gca,'FontSize',13);
set(gca,'FontWeight','bold','FontSize',10)
%legend('\delta  \beta=8\pi/16','\delta  \beta=10\pi/16','\delta  \beta=13\pi/16','\delta  \beta=14\pi/16','\delta  \beta=15\pi/16','location','southwest')
legend('\delta  r_a=0.30,  r=0.20','\delta  r_a=0.30,  r=0.25','\delta  r_a=0.30,  r=0.30','\delta  r_a=0.30,  r=0.35','\delta  r_a=0.30,  r=0.40','location','southwest')
xlabel('\alpha','fontsize',15);
ylabel('\upsilon_{FD}/\upsilon','fontsize',15);

title('TE-C,  \beta=14\pi/6','FontWeight','bold','FontSize',14)

figure;
for ir=1:5
    r=r_stencil(ir);
    rt=0.3;
    %c=fdcoeff_time_space_angles(M,2,v,dt,dx);
    c=fdcoeff_time_space_angles_r(M,0,rt)
    [m n]=size(beta);
    [nn m]=size(theta);
    delta=zeros(n,m);
    for i=1:n    %beta
        bb=beta(i);
        for j=1:m     %theta
            tt=theta(j);
            temp=c(1);
             for k=1:M
               temp=temp+2*c(1+k)*(   cos(k*bb*cos(tt)) + cos(k*bb*sin(tt))          );
             end
            delta(i,j)=1/r/bb*acos( 1+  0.5*r^2*temp              );
        end
    end
    hold on;
    plot(beta,delta(:,1),'LineWidth',1.6)     %  0/16

end

grid on;
set(gca,'linewidth',1.2);
box on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
%axis([0,pi,0.89,1.025])
axis([0,pi,0.91,1.045])
set(gca,'FontSize',13);
set(gca,'FontWeight','bold','FontSize',10)
legend('\delta  r_a=0.30,  r=0.20','\delta  r_a=0.30,  r=0.25','\delta  r_a=0.30,  r=0.30','\delta  r_a=0.30,  r=0.35','\delta  r_a=0.30,  r=0.40','location','southwest')
xlabel('\beta','fontsize',15);
ylabel('\upsilon_{FD}/\upsilon','fontsize',15);
title('TE-C,  \alpha=0','FontWeight','bold','FontSize',14)
