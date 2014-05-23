function [ t, x, z ]=randomModel...
                 (T,Ts,thV,phV,ps0,psV,L,Lft,m,kU,U,noise,linVar,gyroBias0)
% model produces system state and observations for a 2axis chaotic pendulum
%   Models state as x = [th, thd, ph, phd, gyroBias1, gyroBias2]
%   Observations,z =[accX gyroY accY gyroX accZ gyroZ mag mag_1Hz* gps] 
%                                                *interpolated to 100Hz
%   noise is the sensor noise to be modellled
% T - length of time to run model, Ts - sampling period, 
% thV - standard deviation of the random angular motion
% phV - standard deviation of the random angular motion
% Lft - pendulum vertical driving force, m - pendulum mass
% kU/kV - pendulum drag for U&V - wind; for calculating horiz. drive force
% noise - noise to add to the sensors
% accVar - standard deviation of the trans motion
% gyroBias0 - init bias in model. if gyroBias=0, no bias is applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M J HOBBY (2013) mhobby1979@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81; %gravitational constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYSTEM STATE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0:Ts:T;
N=length(t);       
Fu=kU*(U^2);
T=sqrt(Fu^2+(Lft-(m*g))^2);
theta=asin(sign(U)*Fu/T);
th=zeros(N,1); ph=zeros(N,1); ps=zeros(N,1);
nsT=1000;
th(1:nsT*floor(N/nsT))=theta+interp(thV*randn(floor(N/nsT),1),nsT);
for i=nsT*floor(N/nsT):N; th(i)=th(nsT*floor(N/nsT)); end;
thd=[0;diff(th)/Ts];
ph(1:nsT*floor(N/nsT))=interp(phV*randn(floor(N/nsT),1),nsT);
for i=nsT*floor(N/nsT):N; ph(i)=ph(nsT*floor(N/nsT)); end;
phd=[0;diff(ph)/Ts];
nsPS=100;
ps(1:nsPS*floor(N/nsPS))=ps0+interp(psV*randn(floor(N/nsPS),1),nsPS);
for i=nsPS*floor(N/nsPS):N; ps(i)=ps(nsPS*floor(N/nsPS)); end;
psd=[0;diff(ps)/Ts];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACCELERATION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if max(size(linVar))~=3
    linVar=linVar(1)*ones(3,1);
end

xd_e3=zeros(N,1);
xd_e3(1:nsT*floor(N/nsT))=interp(linVar(1)*randn(floor(N/nsT),1),nsT);
for i=nsT*floor(N/nsT):N; xd_e3(i)=xd_e3(nsT*floor(N/nsT)); end;
xdd_e3=[0;diff(xd_e3)/Ts];

yd_e3=zeros(N,1);
yd_e3(1:nsT*floor(N/nsT))=interp(linVar(2)*randn(floor(N/nsT),1),nsT);
for i=nsT*floor(N/nsT):N; yd_e3(i)=yd_e3(nsT*floor(N/nsT)); end;
ydd_e3=[0;diff(yd_e3)/Ts];

zd_e3=zeros(N,1);
zd_e3(1:nsT*floor(N/nsT))=interp(linVar(3)*randn(floor(N/nsT),1),nsT);
for i=nsT*floor(N/nsT):N; zd_e3(i)=zd_e3(nsT*floor(N/nsT)); end;
zdd_e3=[0;diff(zd_e3)/Ts];

xdd_g=-g*sin(th);
ydd_g=g*cos(th).*sin(ph);
zdd_g=g*cos(th).*cos(ph);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GYRO DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%body frame angular rate data
phd_e3=phd-psd.*sin(th);
thd_e3=thd.*cos(ph)+psd.*cos(th).*sin(ph);
psd_e3=-thd.*sin(ph)+psd.*cos(th).*cos(ph);

%gyro drift bias'
bias1=gyroBias0+(Ts*cumtrapz(deg2rad(.00055)*randn(N,1))); 
bias2=gyroBias0+(Ts*cumtrapz(deg2rad(.00055)*randn(N,1))); 
bias3=gyroBias0+(Ts*cumtrapz(deg2rad(.00055)*randn(N,1))); 
if gyroBias0>0
    gyroBias1=bias1;
    gyroBias2=bias2;
    gyroBias3=bias3;
else
    gyroBias1=zeros(N,1);
    gyroBias2=zeros(N,1);
    gyroBias3=zeros(N,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAGNETOMETER DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=[18.25*ones(1,N);0.71*ones(1,N);-45.74*ones(1,N)];
[Bx_e1, By_e1, Bz_e1]=yaw(B, ps, N);
[Bx_e2, By_e2, Bz_e2]=pitch([Bx_e1;By_e1;Bz_e1], th, N);
[Bx_e3, By_e3, Bz_e3]=roll([Bx_e2;By_e2;Bz_e2], ph, N);
%convert to 1Hz data
tBxo=Bx_e3(1:100:N); tBxo=tBxo'+(noise(7)*randn(length(tBxo),1));
Bx_e3I=interp1(1:100:N, tBxo, 1:N);
tByo=By_e3(1:100:N); tByo=tByo'+(noise(8)*randn(length(tByo),1));
By_e3I=interp1(1:100:N, tByo, 1:N);
tBzo=Bz_e3(1:100:N); tBzo=tBzo'+(noise(9)*randn(length(tBzo),1));
Bz_e3I=interp1(1:100:N, tBzo, 1:N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xd_e2, yd_e2, zd_e2]=roll([xd_e3';yd_e3';zd_e3'], -ph, N);
[xd_e1, yd_e1, zd_e1]=pitch([xd_e2;yd_e2;zd_e2], -th, N);
[xd, yd, zd]=yaw([xd_e1;yd_e1;zd_e1], -ps, N);
xdd=[0;diff(xd')/Ts];
ydd=[0;diff(yd')/Ts];
zdd=[0;diff(zd')/Ts];

txd=xd(1:100:N); GPSxd=txd+(noise(10)*randn(1,length(txd)));
GPSxdd=interp1(51:100:N, diff(GPSxd'), 1:N, 'linear', 'extrap');
GPSxd=interp1(1:100:N, GPSxd', 1:N, 'linear', 'extrap');

tyd=yd(1:100:N); GPSyd=tyd+(noise(11)*randn(1,length(tyd)));
GPSydd=interp1(51:100:N, diff(GPSyd'), 1:N, 'linear', 'extrap');
GPSyd=interp1(1:100:N, GPSyd', 1:N, 'linear', 'extrap');

tzd=zd(1:100:N); GPSzd=tzd+(noise(12)*randn(1,length(tzd)));
GPSzdd=interp1(51:100:N, diff(GPSzd'), 1:N, 'linear', 'extrap');
GPSzd=interp1(1:100:N, GPSzd', 1:N, 'linear', 'extrap');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WIND DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(N,1); v=zeros(N,1); w=zeros(N,1);
nsU=500;
U_=5; V_=5;
Uprime=1; Vprime=1; Wprime=1;

u(1:nsU*floor(N/nsU))=U_+interp(Uprime*randn(floor(N/nsU),1),nsU);
for i=nsU*floor(N/nsU):N; u(i)=u(nsU*floor(N/nsU)); end;

v(1:nsU*floor(N/nsU))=V_+interp(Vprime*randn(floor(N/nsU),1),nsU);
for i=nsU*floor(N/nsU):N; v(i)=v(nsU*floor(N/nsU)); end;

w(1:nsU*floor(N/nsU))=interp(Wprime*randn(floor(N/nsU),1),nsU);
for i=nsU*floor(N/nsU):N; w(i)=w(nsU*floor(N/nsU)); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLATE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z(1,:)=xdd_g+xdd_e3+(noise(1)*randn(N,1));                 %x acc
z(2,:)=thd_e3+(noise(2)*randn(N,1))+gyroBias1;             %pitchd
z(3,:)=ydd_g+ydd_e3+(noise(3)*randn(N,1));                 %y acc
z(4,:)=phd_e3+(noise(4)*randn(N,1))+gyroBias2;             %rolld
z(5,:)=zdd_g+zdd_e3+(noise(5)*randn(N,1));                 %z acc
z(6,:)=psd_e3+(noise(6)*randn(N,1))+gyroBias3;
z(7,:)=Bx_e3'+(noise(7)*randn(N,1));                        %Bx
z(8,:)=By_e3'+(noise(8)*randn(N,1));                        %By
z(9,:)=Bz_e3'+(noise(9)*randn(N,1));                        %Bz
z(10,:)=Bx_e3I;                                            %Bx1Hz
z(11,:)=By_e3I;                                            %By1Hz
z(12,:)=Bz_e3I;                                            %Bz1Hz
z(13,:)=GPSxd;
z(14,:)=GPSxdd;
z(15,:)=GPSyd;
z(16,:)=GPSydd;
z(17,:)=GPSzd;
z(18,:)=GPSzdd;
z(19,:)=u-xd_e3;
z(20,:)=v-yd_e3;
z(21,:)=w-zd_e3;
x=[th,thd,ph,phd,ps,psd,gyroBias1,gyroBias2,gyroBias3,...
    xdd,ydd,zdd,xdd_e3,ydd_e3,zdd_e3,u,v,w]; %state
end

