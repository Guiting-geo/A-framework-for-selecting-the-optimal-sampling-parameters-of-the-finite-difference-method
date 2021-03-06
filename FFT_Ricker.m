clc;
clear all;

nt=1024;
fm=35;
dt=0.001;
for it=0:nt-1
    a=pi*fm*(it*dt-1.0/fm);
    a=a*a;
    Ricker(it+1)=(1.0-2.0*a)*exp(-a);
end

fs=1/dt; %采样频率
Ndata=nt; %数据长度
N=512; %FFT的数据长度
n=0:Ndata-1;t=n/fs;   %数据对应的时间序列
y=fft(Ricker,N);   %信号的Fourier变换
mag=abs(y);    %求取振幅
f=(0:N-1)*fs/N; %真实频率
figure;
plot(f(1:80),mag(1:80)/N,'k','LineWidth',1.5); %绘出Nyquist频率之前的振幅

set(gca,'linewidth',1.2);
box on
set(gca,'GridLineStyle','--','GridColor','k','GridAlpha',1);
set(gca,'FontSize',13);
set(gca,'FontWeight','bold','FontSize',10)
%legend('Fixed a','Variable a','Reference','','Reference')
xlabel('f(Hz)','fontsize',15);
ylabel('Amplitude','fontsize',15);
set(gca,'yticklabel','');

