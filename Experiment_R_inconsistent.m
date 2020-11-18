clc;
clear all;

nt=2001;

nx=400;
nz=400;

dt=0.001;
dx=6;
fm=35;

Ricker=zeros(1,nt);
v=ones(nz,nx)*1500;
u0=zeros(nz,nx);
u1=zeros(nz,nx);
u2=zeros(nz,nx);

u0s=zeros(nz,nx);
u1s=zeros(nz,nx);
u2s=zeros(nz,nx);

u0b=zeros(nz,nx);
u1b=zeros(nz,nx);
u2b=zeros(nz,nx);

ckp_0=zeros(20,nx,nz);
ckp_s=ckp_0;
ckp_b=ckp_0;

r=v(1,1)*dt/dx;
for it=0:nt-1
    a=pi*fm*(it*dt-1.0/fm);
    a=a*a;
    Ricker(it+1)=(1.0-2.0*a)*exp(-a);
    
end

M=8;

c=fdcoeff_time_space_angles_r(M,0,r);
cs=fdcoeff_time_space_angles_r(M,0,r*0.6);
cb=fdcoeff_time_space_angles_r(M,0,r*1.4);

MN=M;
index=0;
for it=1:nt
    [1 it]
    u1(nz/2,nx/2)=u1(nz/2,nx/2)+Ricker(it);
    u1s(nz/2,nx/2)=u1s(nz/2,nx/2)+Ricker(it);
    u1b(nz/2,nx/2)=u1b(nz/2,nx/2)+Ricker(it);
    
    for i=MN+1:nz-MN-1
        for j=MN+1:nx-MN-1
            
            diff1=c(1)*u1(i,j);
            diff1s=cs(1)*u1s(i,j);
            diff1b=cb(1)*u1b(i,j);
            for m=1:M
                diff1=diff1+c(m+1)*( u1(i,j+m)+u1(i,j-m)+u1(i+m,j)+u1(i-m,j) );
                diff1s=diff1s+cs(m+1)*( u1s(i,j+m)+u1s(i,j-m)+u1s(i+m,j)+u1s(i-m,j) );
                diff1b=diff1b+cb(m+1)*( u1b(i,j+m)+u1b(i,j-m)+u1b(i+m,j)+u1b(i-m,j) );
            end
            
            
            u2(i,j)=2*u1(i,j)-u0(i,j)+v(i,j)^2*dt^2/dx^2*( diff1 );
            u2s(i,j)=2*u1s(i,j)-u0s(i,j)+v(i,j)^2*dt^2/dx^2*( diff1s );
            u2b(i,j)=2*u1b(i,j)-u0b(i,j)+v(i,j)^2*dt^2/dx^2*( diff1b );
            
        end
    end
    u0=u1;
    u1=u2;
    
    u0s=u1s;
    u1s=u2s;
    
    u0b=u1b;
    u1b=u2b;
    if mod(it,20)==0
            index=index+1;
            ckp_0(index,:,:)=u1;
            ckp_s(index,:,:)=u1s;
            ckp_b(index,:,:)=u1b;
    end
    if it==100
        sp1=u1;
        sp1s=u1s;
        sp1b=u1b;
    end
    if it==300
        sp2=u1;
        sp2s=u1s;
        sp12b=u1b;
    end
    if it==600
        sp3=u1;
        sp3s=u1s;
        sp3b=u1b;
    end
    if it==700
        sp4=u1;
        sp4s=u1s;
        sp4b=u1b;
    end
    if it==1000
        sp5=u1;
        sp5s=u1s;
        sp5b=u1b;
    end
    if it==1200
        sp6=u1;
        sp6s=u1s;
        sp6b=u1b;
    end
    if it==1400
        sp7=u1;
        sp7s=u1s;
        sp7b=u1b;
    end
    if it==1600
        sp8=u1;
        sp8s=u1s;
        sp8b=u1b;
    end
    if it==1800
        sp9=u1;
        sp9s=u1s;
        sp9b=u1b;
    end
    re1(it)=u1(256,256);
    re2(it)=u1(374,207);
    re3(it)=u1(347,165);
    
    re1s(it)=u1s(256,256);
    re2s(it)=u1s(374,207);
    re3s(it)=u1s(347,165);
    
    re1b(it)=u1b(256,256);
    re2b(it)=u1b(374,207);
    re3b(it)=u1b(347,165);
end

figure;
imagesc(sp4);
colormap(gray)
caxis([-0.05 0.05]);
colorbar

fid=fopen('snap_ts_m.bin','wb');
fwrite(fid,sp7,'float32');

fid=fopen('snap_ts_s.bin','wb');
fwrite(fid,sp7s,'float32');

fid=fopen('snap_ts_b.bin','wb');
fwrite(fid,sp7b,'float32');

%{
figure;
imagesc(sp1);
colormap(gray)
caxis([-0.01 0.01]);
colorbar

figure;
imagesc(u1);
colormap(gray)
caxis([-0.1 0.1]);
colorbar

figure; plot(sp1(:,256))
fid=fopen('snap_ts_m.bin','wb');
fwrite(fid,u1,'float32');

fid=fopen('snap_ts_s.bin','wb');
fwrite(fid,u1s,'float32');

fid=fopen('snap_ts_b.bin','wb');
fwrite(fid,u1b,'float32');

%}


