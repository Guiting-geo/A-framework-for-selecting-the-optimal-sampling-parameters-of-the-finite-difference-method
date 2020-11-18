function  a=fdcoeff_time_space_angles_r(M,tt,r)
a=zeros(M+1,1);
A=zeros(M,M);
B=zeros(M,1);

for j=1:M
     for m=1:M
            A(j,m)=m^(2*j)*( (cos(      tt*pi/16     ))^(2*j)+ (sin(  tt*pi/16     ))^(2*j)  );
     
     end
    B(j)=r^(2*j-2);
end

size(A);
size(B);

c=A\B;

a(1)=-4*sum(c);
a(2:M+1)=c;

end
