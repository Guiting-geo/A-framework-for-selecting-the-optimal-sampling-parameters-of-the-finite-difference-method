function  a=fdcoeff_time_space(M,v,dt,dx)
r=v*dt/dx;
a=zeros(M+1,1);

for m=1:M

        a(m+1)=(-1)^(m+1)/m^2;
        temp1=1;
        for i=1:M
                if i ~= m
                    temp1=temp1*abs(  (i^2-r^2)/(i^2-m^2)     );
                end
        end
        
        a(m+1)=a(m+1)*temp1; 

end

    a(1)=-4*sum(a(2:end));

end
