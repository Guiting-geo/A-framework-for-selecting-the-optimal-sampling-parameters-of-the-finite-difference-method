function r_tste=stability_tste_r_M_plot(M)

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


for i=1:1000
    r(i)=0.001*i;
    a=fdcoeff_time_space_angles_r(M,0,r(i));
    
    temp=0;
    for m=1:M
        temp=temp+a(m+1)*(  (-1)^(m-1)  +1   );
    end
    s_tste(i)=1/sqrt(temp);
end

s_max=zeros(1,12);

for im=1:12
    
    for i=1:1000
        r(i)=0.001*i;
        a=fdcoeff_time_space_angles_r(im,0,r(i));
        
        temp=0;
        for m=1:im
            temp=temp+a(m+1)*(  (-1)^(m-1)  +1   );
        end
        s_m(i)=1/sqrt(temp);
        if s_m(i)<r(i)
            s_max(im)=s_m(i);
            break;
        end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(r(1:i_tse),s_tste(1:i_tse),'m','LineWidth',3.5); hold on   %  admissible s(r)




plot(r,s_tste,'--g','LineWidth',2.0); hold on         % s(r)

R=0.01:0.01:1;
plot(R,R,'b','LineWidth',1.6); hold on          %%%%  r=r


%%%%%%%%%%%%%%%%%%%%%%%% A, B, C points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(r(i_tse),s_tste(i_tse),'or','LineWidth',1.6); hold on         %%% A point

text(r(i_tse),s_tste(i_tse)-0.01,'A','FontSize',13);                 % A


XRmin=ones(1,100)*s_tste(i_tse);
YRmin=(1:100)*0.01;
size(XRmin);
plot(XRmin(1:5:end),YRmin(1:5:end),'--','Color',[1 0.5 0],'LineWidth',1.6); hold on         %%% Rmin
text(s_tste(i_tse)-0.03,0.1,'R_{max}','FontSize',13);

%%%%%%%%%%%%%%%%%%%  variation of velocity  %%%%%%%%%%%%%%%%%%%

axis([0.,1,0,1])
set(gca,'linewidth',1.2);
box on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
set(gca,'FontSize',13);

set(gca,'FontWeight','bold','FontSize',12)
legend('Admissible s(r_a)','s(r_a)','s(r_a)=r_a','location','northwest')

xlabel('r_a','fontsize',15);
ylabel('s(r_a)','fontsize',15);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  M  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MM=1:12;
figure;
plot(MM,s_max,'-or','LineWidth',1.6); hold on         %%% A point

set(gca,'linewidth',1.2);
box on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
set(gca,'FontSize',13);

set(gca,'FontWeight','bold','FontSize',12)
legend('Point A','location','northeast')

xlabel('M','fontsize',15);
ylabel('R_{max}','fontsize',15);



end




















