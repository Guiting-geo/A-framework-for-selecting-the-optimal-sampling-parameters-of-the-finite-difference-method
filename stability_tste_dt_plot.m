function dt=stability_tste_dt_plot(M,Vmin,Vmax, dx)

for i=1:1000
    r(i)=0.001*i;
    a=fdcoeff_time_space_angles_r(M,0,r(i));
    temp=0;
    for m=1:M
        temp=temp+a(m+1)*(  (-1)^(m-1)  +1   );
    end
    s_tste(i)=1/sqrt(temp);
    if s_tste(i)<r(i)
        i_tse=i;
        break;
    end
end
i_tse=i_tse-1;
Ra=(i_tse)*0.001;

tmax=Ra*dx/Vmax;


dtt=0.00005:0.000001:tmax;
[m n]=size(dtt);
i_op=0;

minic=18000000;
for i=1:n
    
    dt=dtt(i);
    Rmin=Vmin*dt/dx;
    Rmax=Vmax*dt/dx;
    a=fdcoeff_time_space_angles_r(M,0,Rmin);
    temp=0;
    for m=1:M
        temp=temp+a(m+1)*(  (-1)^(m-1)  +1   );
    end
    s_ts(i)=1/sqrt(temp);
    if abs(s_ts(i)-Rmax)<=minic
        i_op=i;
        minic=abs(s_ts(i)-Rmax);
    end
    
end

[i_tse i_op];
Rmax=Vmax*dtt(i_op)/dx;
Rmin=Vmin*dtt(i_op)/dx;
if i_op>0
    
    [dtt(i_op)  Rmin Rmax s_ts(i_op)]
    
    figure;
    irmin=round(dtt(i_op)*Vmin/dx/0.001);
    irmax=floor(Rmax/0.001);
    [irmin irmax];
    plot(r(irmin:1:irmax),s_tste(irmin:1:irmax),'m','LineWidth',3.5); hold on   %  optimal s(r)
    
    
    plot(r,s_tste,'--g','LineWidth',2.0); hold on         % s(r)
    
    R=0.01:0.01:1;
    plot(R,R,'b','LineWidth',1.6); hold on          %%%%  r=r
    
    [m n]=size(r);
    rmax=ones(size(r))*Rmax;
    plot(r,rmax,'r','LineWidth',1.6); hold on          %%%%  rmax
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% A, B, C points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plot(r(end-1),s_tste(end-1),'or','LineWidth',1.6); hold on         %%% A point
    
    plot(s_ts(i_op),s_ts(i_op),'or','LineWidth',1.6); hold on         %%% B point
    
    plot(dtt(i_op)*Vmin/dx,s_ts(i_op),'or','LineWidth',1.6); hold on         %%  C  point
    
    
    text(r(end),s_tste(end)-0.01,'A','FontSize',13);                 % A
    
    text(Rmax,Rmax-0.02,'B','FontSize',13);                          % B
    
    text(dtt(i_op)*Vmin/dx,s_ts(i_op)-0.02,'C','FontSize',13);       % C
    
       
    XRmin=ones(1,100)*dtt(i_op)*Vmin/dx;
    YRmin=(1:100)*0.01;
    size(XRmin);
    plot(XRmin(1:5:end),YRmin(1:5:end),'--','Color',[1 0.5 0],'LineWidth',1.6); hold on         %%% Rmin
    text(dtt(i_op)*Vmin/dx-0.03,0.32,'R_{min}','FontSize',13);
    
    
    XRmax=ones(1,100)*Rmax;
    YRmax=(1:100)*0.01;
    plot(XRmax(1:5:end),YRmax(1:5:end),'--','Color',[1 0.5 0],'LineWidth',1.6); hold on         %%% Rmax
    text(Rmax-0.03,0.32,'R_{max}','FontSize',13);
    
    
    axis([0.,0.6,0.3,0.7]);
    set(gca,'linewidth',1.2);
    box on
    set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
    set(gca,'FontSize',13);
    
    set(gca,'FontWeight','bold','FontSize',9)
    legend('Optimal s(r_a)','s(r_a)','s(r_a)=r_a','s(r_a)=r_{max}','location','northwest')
    
    xlabel('r_a','fontsize',15);
    ylabel('s(r_a)','fontsize',15);
    
end

dt=dtt(i_op);
end
