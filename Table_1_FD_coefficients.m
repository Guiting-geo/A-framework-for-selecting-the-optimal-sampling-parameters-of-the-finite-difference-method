clc;
clear all;

format long

dx=6.;
fm=35;
M=8;

Vmin=1500;
Vm=2500;
Vmax=5500;

dt=stability_tste_dt_plot(M,Vmin,Vmax, dx)

r1=Vmin*dt/dx;
r2=Vm*dt/dx;
r3=Vmax*dt/dx;
a1=fdcoeff_time_space_angles_r(M,0,r1);
a2=fdcoeff_time_space_angles_r(M,0,r2);
a3=fdcoeff_time_space_angles_r(M,0,r3);
