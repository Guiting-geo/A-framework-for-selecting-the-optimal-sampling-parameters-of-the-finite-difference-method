function r_tste=stability_tste_ABC_plot(M,Vmin,Vmax, dt, dx)

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

Rmax=Vmax*dt/dx;
i_op=0;
if Rmax>s_tste(1)
    for i=1:i_tse
        if abs(s_tste(i)-Rmax) <=0.001
            i_op=i;
            break;
        end
    end
end

irmin=round(dt*Vmin/dx/0.001);
irmax=floor(Rmax/0.001);
[irmin*0.001 irmax*0.001 i_op*0.001]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
if i_op>0
    if r(i_op)<dt*Vmin/dx
        plot(r(irmin:1:irmax),s_tste(irmin:1:irmax),'m','LineWidth',3.5); hold on   %  optimal s(r)
    else
        plot(r(i_op:1:irmax),s_tste(i_op:1:irmax),'m','LineWidth',3.5); hold on   %  optimal s(r)
    end
else
    plot(r(irmin:1:irmax),s_tste(irmin:1:irmax),'m','LineWidth',3.5); hold on   %  optimal s(r)
end

plot(r,s_tste,'--g','LineWidth',2.0); hold on         % s(r)

R=0.01:0.01:1;
plot(R,R,'b','LineWidth',1.6); hold on          %%%%  r=r

[m n]=size(r);
rmax=ones(size(r))*Rmax;
plot(r,rmax,'r','LineWidth',1.6); hold on          %%%%  rmax

%%%%%%%%%%%%%%%%%%%%%%%% A, B, C points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(r(end-1),s_tste(end-1),'or','LineWidth',1.6); hold on         %%% A point

plot(Rmax,Rmax,'or','LineWidth',1.6); hold on         %%% B point

text(r(end),s_tste(end)-0.01,'A','FontSize',13);                 % A

text(Rmax,Rmax-0.02,'B','FontSize',13);                          % B

if Rmax>s_tste(1) && i_op>0
    
    plot(r(i_op),s_tste(i_op),'or','LineWidth',1.6); hold on         %%% B point
    text(r(i_op),s_tste(i_op)-0.02,'C','FontSize',13);       % C
end

XRmin=ones(1,100)*dt*Vmin/dx;
YRmin=(1:100)*0.01;
size(XRmin);
plot(XRmin(1:5:end),YRmin(1:5:end),'--','Color',[1 0.5 0],'LineWidth',1.6); hold on         %%% Rmin
text(dt*Vmin/dx-0.03,0.32,'R_{min}','FontSize',13);

XRmax=ones(1,100)*Rmax;
YRmax=(1:100)*0.01;
plot(XRmax(1:5:end),YRmax(1:5:end),'--','Color',[1 0.5 0],'LineWidth',1.6); hold on         %%% Rmax
text(Rmax-0.03,0.32,'R_{max}','FontSize',13);

axis([0.,0.6,0.3,0.7])
set(gca,'linewidth',1.2);
box on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
set(gca,'FontSize',13);

set(gca,'FontWeight','bold','FontSize',9)
legend('Admissible s(r_a)','s(r_a)','s(r_a)=r_a','s(r_a)=R_{max}','location','northwest')

xlabel('r_a','fontsize',15);
ylabel('s(r_a)','fontsize',15);

end




















