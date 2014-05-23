%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL MOTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model instance for 1D motion. 
% estimates of : pitch pitchRate and gyrobias.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
addpath '../general/'
figRoot=...
    'Z:/instrument/Projects/Tethersonde/phd/chapters/chapter6_fig/v4/6_3/';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL CONSTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;
fs=100;
Ts=1/fs; % sampling period
L=100; %pendulum length
m=3.5; %pendulum weight
T=120.01; %model duration (s)
Lft=250; %inverted pendulum lift
kU=.8; %drag coefficient
U=10; %wind speed
B=[18.25,0.71,-45.74];
noise=[.0075,.0008,.0075,.0008,.0075,.0008,.038,.038,.038,.1,.1,.1];%obs noise 
N=(T/Ts)-1; % number of samples expected in model
gyroBiasBar=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODELS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM WITH TRANS MOTION IN X/Y/Z:
[t1,x1,z1]=randomModel(T,Ts,.07,.07,pi/6,.25,L,Lft,m,kU,U,noise,[5 5 1],gyroBiasBar);
    th_=asin(mean(z1(1,:)/-g));                                   
    ph_=asin(mean(z1(3,:))/(g*cos(th_)));                                                    
    x0=[th_ z1(2,1) ph_ z1(4,1) 0 0 z1(13,1) 0 z1(15,1) 0 z1(17,1) 0 0 0 0 0 0 0 0];
    v=[.0075;.0008;.0075;.0008;.0008;.1;.1;.1;.038;.038;.038;.0075;.01;.01;.01];
    w=[10^-6;.1;10^-6;.1;10^-6;.1;10^-6;1;10^-6;1;10^-6;1;.0001;.0001;.0001;1;1;1;1];
    zEKF=z1([1:4,6,13,15,17,7:9,5,19:21],:);
    xhatEKF1=ekf(zEKF,Ts,v,w,x0,B);%100Hz mag obs + accZ + sonicW
	
% PENDULUM:
[t2,x2,z2]=balloonModel(T,Ts,.07,.07,.01,pi/6,.25,L,Lft,m,kU,U,noise,gyroBiasBar);                     
    th_=asin(mean(z2(1,:)/-g));                                   
    ph_=asin(mean(z2(3,:))/(g*cos(th_)));                                                    
    x0=[th_ z2(2,1) ph_ z2(4,1) 0 0 z2(13,1) 0 z2(15,1) 0 z2(17,1) 0 0 0 0 0 0 0 0];
    v=[.0075;.0008;.0075;.0008;.0008;.1;.1;.1;.038;.038;.038;.0075;.01;.01;.01];
    w=[10^-6;.1;10^-6;.1;10^-6;.1;10^-6;1;10^-6;1;10^-6;1;.0001;.0001;.0001;1;1;1;1];
    zEKF=z2([1:4,6,13,15,17,7:9,5,19:21],:);
    xhatEKF2=ekf(zEKF,Ts,v,w,x0,B);%100Hz mag obs + accZ + sonicW
	
% COMBO PENDULUM & RANDOM MOTION: 
[t3,x3,z3]=comboModel(T,Ts,.07,.07,.07,.01,.07,pi/6,.25,[5 5 1],...
                   L,Lft,m,kU,U,noise,gyroBiasBar);   
    th_=asin(mean(z3(1,:)/-g));                                   
    ph_=asin(mean(z3(3,:))/(g*cos(th_)));  
    x0=[th_ z3(2,1) ph_ z3(4,1) 0 0 z3(13,1) 0 z3(15,1) 0 z3(17,1) 0 0 0 0 0 0 0 0];
    v=[.0075;.0008;.0075;.0008;.0008;.1;.1;.1;.038;.038;.038;.0075;.01;.01;.01];
    w=[10^-6;.1;10^-6;.1;10^-6;.1;10^-6;1;10^-6;1;10^-6;1;.0001;.0001;.0001;1;1;1;1];
    zEKF=z3([1:4,6,13,15,17,7:9,5,19:21],:);
    xhatEKF3=ekf(zEKF,Ts,v,w,x0,B);%100Hz mag obs + accZ + sonicW    