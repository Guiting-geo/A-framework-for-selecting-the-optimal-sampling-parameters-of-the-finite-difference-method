clc;
clear all;

nt=1800;
nx=400;
nz=400;
nb=40;

nnx=nx+2*nb;
nnz=nz+2*nb;

dx=6.;
fm=35;
M=12;

Vmin=1500;
Vm=2500;
Vmax=5500;

dt=stability_tste_dt_plot(M,Vmin,Vmax, dx);

nt=floor(1./dt);

Ricker=zeros(1,nt);

v=ones(nz,nx)*Vmin;
v(100:end,:)=Vm;
v(200:end,:)=Vmax;
vv=ones(nnz,nnx)*3500;
vv(nb+1:nb+nz,nb+1:nb+nx)=v(:,:);
for i=1:nb
    vv(i,:)=vv(nb+1,:);
    vv(nb+nz+i,:)=vv(nb+nz-1,:);
    vv(:,i)=vv(:,nb+1);
    vv(:,nb+nx+i)=vv(:,nb+nx-1);
end

for it=0:nt-1
    a=pi*fm*(it*dt-1.0/fm);
    a=a*a;
    Ricker(it+1)=(1.0-2.0*a)*exp(-a);
end

abc_f=zeros(1,nb);

for ib=0:nb-1
    tmp=0.025*(nb-1-ib);
    abc_f(ib+1)=exp(-tmp*tmp);
end

record=zeros(nt,nx);
recordv=zeros(nt,nx);

u0=zeros(nnz,nnx);
u1=zeros(nnz,nnx);
u2=zeros(nnz,nnx);

uv0=zeros(nnz,nnx);
uv1=zeros(nnz,nnx);
uv2=zeros(nnz,nnx);

nck=40;
ckpoints=zeros(nck,nz,nx);
ckpointsv=zeros(nck,nz,nx);
index=0;

cf=zeros(nnz,nnx,M+1);
r=mean(mean(v))*dt/dx;
c=fdcoeff_time_space_angles_r(M,0,r);

for i=1:nnz
    for j=1:nnx
        rr=vv(i,j)*dt/dx;
        cf(i,j,:)=fdcoeff_time_space_angles_r(M,0,rr);
    end
end

for it=1:nt
    [nt it]
    u1(nb+2,nb+200)=u1(nb+2,nb+200)+Ricker(it);
    uv1(nb+2,nb+200)=uv1(nb+2,nb+200)+Ricker(it);
    for i=M+1:nnz-M-1
        for j=M+1:nnx-M-1
            
            diff1=c(1)*u1(i,j);
            diffv1=cf(i,j,1)*uv1(i,j);
            
            for m=1:M
                diffv1=diffv1+cf(i,j,m+1)*( uv1(i,j+m)+uv1(i,j-m)+uv1(i+m,j)+uv1(i-m,j) );
                diff1=diff1+c(m+1)*( u1(i,j+m)+u1(i,j-m)+u1(i+m,j)+u1(i-m,j) );
            end
            uv2(i,j)=2*uv1(i,j)-uv0(i,j)+vv(i,j)^2*dt^2/dx^2*( diffv1);
            u2(i,j)=2*u1(i,j)-u0(i,j)+vv(i,j)^2*dt^2/dx^2*( diff1);
        end
    end
    record1(it) = u1(nb+50,nb+200);
    recordv1(it) =uv1(nb+50,nb+200);
    
    record2(it) = u1(nb+100,nb+200);
    recordv2(it) =uv1(nb+100,nb+200);
    
    record3(it) = u1(nb+150,nb+200);
    recordv3(it) =uv1(nb+150,nb+200);
    
    record4(it) = u1(nb+200,nb+200);
    recordv4(it) =uv1(nb+200,nb+200);
    
    for i=1:nnx
        u0(1:nb,i)=u0(1:nb,i).*abc_f';
        u0(end:-1:end-nb+1,i)=u0(end:-1:end-nb+1,i).*abc_f';
        u1(1:nb,i)=u1(1:nb,i).*abc_f';
        u1(end:-1:end-nb+1,i)=u1(end:-1:end-nb+1,i).*abc_f';
        u2(1:nb,i)=u2(1:nb,i).*abc_f';
        u2(end:-1:end-nb+1,i)=u2(end:-1:end-nb+1,i).*abc_f';
        
        uv0(1:nb,i)=uv0(1:nb,i).*abc_f';
        uv0(end:-1:end-nb+1,i)=uv0(end:-1:end-nb+1,i).*abc_f';
        uv1(1:nb,i)=uv1(1:nb,i).*abc_f';
        uv1(end:-1:end-nb+1,i)=uv1(end:-1:end-nb+1,i).*abc_f';
        uv2(1:nb,i)=uv2(1:nb,i).*abc_f';
        uv2(end:-1:end-nb+1,i)=uv2(end:-1:end-nb+1,i).*abc_f';
    end
    for i=1:nnz
        u0(i,1:nb)=u0(i,1:nb).*abc_f;
        u0(i,end:-1:end-nb+1)=u0(i,end:-1:end-nb+1).*abc_f;
        u1(i,1:nb)=u1(i,1:nb).*abc_f;
        u1(i,end:-1:end-nb+1)=u1(i,end:-1:end-nb+1).*abc_f;
        u2(i,1:nb)=u2(i,1:nb).*abc_f;
        u2(i,end:-1:end-nb+1)=u2(i,end:-1:end-nb+1).*abc_f;
        
        
        uv0(i,1:nb)=uv0(i,1:nb).*abc_f;
        uv0(i,end:-1:end-nb+1)=uv0(i,end:-1:end-nb+1).*abc_f;
        uv1(i,1:nb)=uv1(i,1:nb).*abc_f;
        uv1(i,end:-1:end-nb+1)=uv1(i,end:-1:end-nb+1).*abc_f;
        uv2(i,1:nb)=uv2(i,1:nb).*abc_f;
        uv2(i,end:-1:end-nb+1)=uv2(i,end:-1:end-nb+1).*abc_f;
        
    end
    if mod(it,40)==0 && it>800 && index <= nck
        index=index+1;
        ckpoints(index,:,:)=u1(nb+1:nb+nz,nb+1:nb+nx);
        ckpointsv(index,:,:)=uv1(nb+1:nb+nz,nb+1:nb+nx);
    end
    
    u0=u1;
    u1=u2;
    
    uv0=uv1;
    uv1=uv2;
end

%{
sp=reshape(ckpoints(15,:,:),nz,nx);
figure;
imagesc(sp);
colormap(gray)
caxis([-0.1 0.1])

%fid=fopen('snap_ts_fixed_v.bin','wb');
%fwrite(fid,sp,'float32');


sp=reshape(ckpointsv(15,:,:),nz,nx);
figure;
imagesc(sp);
colormap(gray)
caxis([-0.1 0.1])

%fid=fopen('snap_ts_variable_v.bin','wb');
%fwrite(fid,sp,'float32');


r1=record3(800:1600)*0.2;
r2=recordv3(800:1600)*0.2;
t=(800:1600)*dt;

figure;
plot(t,r1); hold on;
plot(t,r2-0.2); hold on;

axis([0.48,0.87,-0.35,0.22])
%}

figure;
imagesc(u1(nb+1:nb+nz,nb+1:nb+nx));
colormap(gray)
caxis([-0.1 0.1])

figure;
imagesc(uv1(nb+1:nb+nz,nb+1:nb+nx));
colormap(gray)
caxis([-0.1 0.1])

fid=fopen('snap_ts_fixed_v_M12.bin','wb');
fwrite(fid,u1(nb+1:nb+nz,nb+1:nb+nx),'float32');

fid=fopen('snap_ts_variable_v_M12.bin','wb');
fwrite(fid,uv1(nb+1:nb+nz,nb+1:nb+nx),'float32');
%fwrite(fid,sp,'float32');
