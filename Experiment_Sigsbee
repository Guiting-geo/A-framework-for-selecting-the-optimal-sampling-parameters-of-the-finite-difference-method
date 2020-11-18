clc;
clear all;

fid=fopen('sigsbee_1201_3201.bin','rb');
sig=fread(fid,[1201 3201],'float32');
model=sig(1:1000,800:end)*1000;
[mm nn]=size(model);

nz=floor(mm*0.4);
nx=floor(nn*0.4);

vt=imresize(model,[nz nx],'bilinear');
Vmin=min(min(vt));
Vmax=max(max(vt));
vvmin=1550;
vvmax=4500;
dsc=Vmax-Vmin;
dsc2=vvmax-vvmin;
v=(vt-Vmin)/dsc*dsc2+vvmin;

fidv=fopen('Sig_model.bin','wb');
fwrite(fidv,v,'float32');

nt=3500;
nb=40;

nnx=nx+2*nb;
nnz=nz+2*nb;

%dt=0.0007;
dx=6.;
fm=40;
Ricker=zeros(1,nt);
M=8;
dt=stability_tste_dt_plot(M,vvmin,vvmax, dx);

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

nck=(nt-1000)/100;
ckpoints=zeros(nck,nz,nx);
ckpointsv=zeros(nck,nz,nx);
index=0;

M=8;

cf=zeros(nnz,nnx,M+1);
r=0.36;
c=fdcoeff_time_space_angles_r(M,0,r);

for i=1:nnz
    for j=1:nnx
        rr=vv(i,j)*dt/dx;
        cf(i,j,:)=fdcoeff_time_space_angles_r(M,0,rr);
    end
end

for it=1:nt
    [2 it]
    u1(nb+2,nb+nx/2)=u1(nb+2,nb+nx/2)+Ricker(it);
    uv1(nb+2,nb+nx/2)=uv1(nb+2,nb+nx/2)+Ricker(it);
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
    record(it,:)=u1(nb+5,nb+1:nx+nb);
    recordv(it,:)=uv1(nb+5,nb+1:nx+nb);
    
    signal1(it)=u0(nb+nz/4,nb+nx/2);
    signal2(it)=u0(nb+nz/2,nb+nx/2);
    signal3(it)=u0(nb+nz*3/4,nb+nx/2);
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
    u0=u1;
    u1=u2;
    
    uv0=uv1;
    uv1=uv2;
    if mod(it,100)==0 && it>1000 && index <= nck
        index=index+1;
        ckpoints(index,:,:)=u1(nb+1:nb+nz,nb+1:nb+nx);
        ckpointsv(index,:,:)=uv1(nb+1:nb+nz,nb+1:nb+nx);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sp=reshape(ckpoints(10,:,:),nz,nx);
figure;
imagesc(sp);
colormap(gray)
caxis([-0.1 0.1])


sp=reshape(ckpointsv(10,:,:),nz,nx);
figure;
imagesc(sp);
colormap(gray)
caxis([-0.1 0.1])


sp=reshape(ckpointsv(10,:,:),nz,nx);
fids1=fopen('Sig_TS_variable_1.bin','wb');
fwrite(fids1,sp,'float32');


sp=reshape(ckpoints(10,:,:),nz,nx);
fids2=fopen('Sig_TS_fixed_1.bin','wb');
fwrite(fids2,sp,'float32');


fidr1=fopen('Sig_TS_variable_record.bin','wb');
fwrite(fidr1,recordv,'float32');
 
fidr2=fopen('Sig_TS_fixed_record.bin','wb');
fwrite(fidr2,record,'float32');



