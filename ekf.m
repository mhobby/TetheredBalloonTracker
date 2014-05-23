function [ xhat  ] = ekf( z, Ts, v, w, x0, B)
%EKF: 3DOF extended kalman filter solution for producing estimates of
%pitch, roll & yaw using accelerometer, gyroscope, GPS, magnetometer
%     & sonic observations
%  xhat - estimate [ th thd ph phd ps psd xd3 xdd3 yd3 ydd3 zd3 zdd3 bth bph bps zcf u v w]
%     z - obs [ accX gyroY accY gyroX gyroZ GPSXd GPSYd GPSZd Bx By Bz accZ sonicU sonicV sonicW]
%    Ts - sample period
%   n_v - vector (1x19) of measurement noise 
%   n_w - vector (1x15) of process noise 
%    x0 - initial estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M J HOBBY (2013) mhobby1979@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=9.81;
N=length(z);

xhat=zeros(19,N);
hxk_pre=zeros(15,1);

G=[ 1  0    0  0    0  0    0  0    0  0    0  0    0  0  0  0  0  0  0; %th
    0  1    0  0    0  0    0  0    0  0    0  0    0  0  0  0  0  0  0; %thd
    0  0    1  0    0  0    0  0    0  0    0  0    0  0  0  0  0  0  0; %ph
    0  0    0  1    0  0    0  0    0  0    0  0    0  0  0  0  0  0  0; %phd
    0  0    0  0    1  0    0  0    0  0    0  0    0  0  0  0  0  0  0; %ps
    0  0    0  0    0  1    0  0    0  0    0  0    0  0  0  0  0  0  0; %psd    
    0  0    0  0    0  0    1  0    0  0    0  0    0  0  0  0  0  0  0; %xd
    0  0    0  0    0  0    0  1    0  0    0  0    0  0  0  0  0  0  0; %xdd
    0  0    0  0    0  0    0  0    1  0    0  0    0  0  0  0  0  0  0; %yd
    0  0    0  0    0  0    0  0    0  1    0  0    0  0  0  0  0  0  0; %ydd
    0  0    0  0    0  0    0  0    0  0    1  0    0  0  0  0  0  0  0; %zd
    0  0    0  0    0  0    0  0    0  0    0  1    0  0  0  0  0  0  0; %zdd
    0  0    0  0    0  0    0  0    0  0    0  0   Ts  0  0  0  0  0  0; %bth
    0  0    0  0    0  0    0  0    0  0    0  0    0 Ts  0  0  0  0  0; %bph
    0  0    0  0    0  0    0  0    0  0    0  0    0  0 Ts  0  0  0  0; %bps
    0  0    0  0    0  0    0  0    0  0    0  0    0  0  0  1  0  0  0; %zcf
    0  0    0  0    0  0    0  0    0  0    0  0    0  0  0  0  1  0  0; %u
    0  0    0  0    0  0    0  0    0  0    0  0    0  0  0  0  0  1  0; %v
    0  0    0  0    0  0    0  0    0  0    0  0    0  0  0  0  0  0  1];%w

w=[w(1);w(2);w(3);w(4);w(5);w(6);w(7);...
   w(8);w(9);w(10);w(11);w(12);w(13);w(14);w(15);w(16);w(17);w(18);w(19)];
Q=diag(diag((G*w)*(G*w)'));                    % process covariance matrix

R=[v(1)^2 0    0    0    0    0    0    0    0     0     0     0     0    0     0   ;
     0  v(2)^2 0    0    0    0    0    0    0     0     0     0     0    0     0   ;
     0    0  v(3)^2 0    0    0    0    0    0     0     0     0     0    0     0   ;
     0    0    0  v(4)^2 0    0    0    0    0     0     0     0     0    0     0   ;
     0    0    0    0  v(5)^2 0    0    0    0     0     0     0     0    0     0   ;
     0    0    0    0    0  v(6)^2 0    0    0     0     0     0     0    0     0   ;
     0    0    0    0    0    0  v(7)^2 0    0     0     0     0     0    0     0   ;
     0    0    0    0    0    0    0  v(8)^2 0     0     0     0     0    0     0   ;
     0    0    0    0    0    0    0    0  v(9)^2  0     0     0     0    0     0   ;
     0    0    0    0    0    0    0    0    0  v(10)^2  0     0     0    0     0   ;
     0    0    0    0    0    0    0    0    0     0  v(11)^2  0     0    0     0   ;
     0    0    0    0    0    0    0    0    0     0     0  v(12)^2  0    0     0   ;
     0    0    0    0    0    0    0    0    0     0     0     0  v(13)^2 0     0   ;
     0    0    0    0    0    0    0    0    0     0     0     0     0  v(14)^2 0   ;
     0    0    0    0    0    0    0    0    0    0     0     0     0     0  v(15)^2];
xhat(:,1)=x0;

% initial estimate covariance: 
P_post(:,:,1)=Q;
for i=2:N
    F = [ 1 Ts    0  0    0  0    0  0   0  0   0  0  0  0  0  0 0 0 0; %th   1
          0  1    0  0    0  0    0  0   0  0   0  0  0  0  0  0 0 0 0; %thd  2
          0  0    1 Ts    0  0    0  0   0  0   0  0  0  0  0  0 0 0 0; %ph   3
          0  0    0  1    0  0    0  0   0  0   0  0  0  0  0  0 0 0 0; %phd  4
          0  0    0  0    1 Ts    0  0   0  0   0  0  0  0  0  0 0 0 0; %ps   5
          0  0    0  0    0  1    0  0   0  0   0  0  0  0  0  0 0 0 0; %psd  6          
          0  0    0  0    0  0    1 Ts   0  0   0  0  0  0  0  0 0 0 0; %xd3  7
          0  0    0  0    0  0    0  1   0  0   0  0  0  0  0  0 0 0 0; %xdd3 8
          0  0    0  0    0  0    0  0   1 Ts   0  0  0  0  0  0 0 0 0; %yd3  9
          0  0    0  0    0  0    0  0   0  1   0  0  0  0  0  0 0 0 0; %ydd3 10
          0  0    0  0    0  0    0  0   0  0   1 Ts  0  0  0  0 0 0 0; %zd3  11
          0  0    0  0    0  0    0  0   0  0   0  1  0  0  0  0 0 0 0; %zdd3 12                   
          0  0    0  0    0  0    0  0   0  0   0  0  1  0  0  0 0 0 0; %bth  13
          0  0    0  0    0  0    0  0   0  0   0  0  0  1  0  0 0 0 0; %bph  14
          0  0    0  0    0  0    0  0   0  0   0  0  0  0  1  0 0 0 0; %bps  15
          0  0    0  0    0  0    0  0   0  0   0  0  0  0  0  1 0 0 0; %zcf  16
          0  0    0  0    0  0    0  0   0  0   0  0  0  0  0  0 1 0 0; %u    17
          0  0    0  0    0  0    0  0   0  0   0  0  0  0  0  0 0 1 0; %v    18
          0  0    0  0    0  0    0  0   0  0   0  0  0  0  0  0 0 0 1];%w    19
    xhat_pre=F*xhat(:,i-1);
    P_pre=F*P_post*F' + Q; 
        
    hxk_pre(1)=-g*sin(xhat_pre(1))+xhat_pre(8);  
    hxk_pre(2)=(xhat_pre(2)*cos(xhat_pre(3)))...
              +(xhat_pre(6)*cos(xhat_pre(1))*sin(xhat_pre(3)))+xhat_pre(13);
    hxk_pre(3)=g*cos(xhat_pre(1))*sin(xhat_pre(3))+xhat_pre(10);
    hxk_pre(4)=xhat_pre(4)-(xhat_pre(6)*sin(xhat_pre(1)))+xhat_pre(14);      
    hxk_pre(5)=-(xhat_pre(2)*sin(xhat_pre(3)))...
              +(xhat_pre(6)*cos(xhat_pre(1))*cos(xhat_pre(3)))...
              +xhat_pre(15);
    hxk_pre(6)=(xhat_pre(7)*cos(xhat_pre(1))*cos(xhat_pre(5)))...
              +(xhat_pre(9)*(-sin(xhat_pre(5))*cos(xhat_pre(3))...
                    +sin(xhat_pre(1))*sin(xhat_pre(3))*cos(xhat_pre(5))) )...
              +(xhat_pre(11)*(sin(xhat_pre(3))*sin(xhat_pre(5))...
                    +sin(xhat_pre(1))*cos(xhat_pre(3))*cos(xhat_pre(5))) );
    hxk_pre(7)=(xhat_pre(7)*sin(xhat_pre(5))*cos(xhat_pre(1)))...
              +(xhat_pre(9)*(cos(xhat_pre(3))*cos(xhat_pre(5))...
                    +sin(xhat_pre(1))*sin(xhat_pre(3))*sin(xhat_pre(5))) )...
              +(xhat_pre(11)*(-sin(xhat_pre(3))*cos(xhat_pre(5))...
                    +sin(xhat_pre(1))*cos(xhat_pre(3))*sin(xhat_pre(5))) );
    hxk_pre(8)=(-xhat_pre(7)*sin(xhat_pre(1)))...
              +(xhat_pre(9)*cos(xhat_pre(1))*sin(xhat_pre(3)))...
              +(xhat_pre(11)*cos(xhat_pre(1))*cos(xhat_pre(3)));
    hxk_pre(9)=B(1)*cos(xhat_pre(5))*cos(xhat_pre(1))...
              +B(2)*sin(xhat_pre(5))*cos(xhat_pre(1))...
              -B(3)*sin(xhat_pre(1));
    hxk_pre(10)=B(1)*(-sin(xhat_pre(5))*cos(xhat_pre(3))...
                    +cos(xhat_pre(5))*sin(xhat_pre(3))*sin(xhat_pre(1)))...
               +B(2)*(cos(xhat_pre(5))*cos(xhat_pre(3))...
                    +sin(xhat_pre(5))*sin(xhat_pre(3))*sin(xhat_pre(1)))...
               +B(3)*cos(xhat_pre(1))*sin(xhat_pre(3));    
    hxk_pre(11)=B(1)*(sin(xhat_pre(5))*sin(xhat_pre(3))...
                    +cos(xhat_pre(5))*cos(xhat_pre(3))*sin(xhat_pre(1)))...
               +B(2)*(-cos(xhat_pre(5))*sin(xhat_pre(3))...
                    +sin(xhat_pre(5))*cos(xhat_pre(3))*sin(xhat_pre(1)))...
               +B(3)*cos(xhat_pre(1))*cos(xhat_pre(3));     
    hxk_pre(12)=g*cos(xhat_pre(1))*cos(xhat_pre(3))...
               +xhat_pre(12)+xhat_pre(16);
    hxk_pre(13)=xhat_pre(17)-xhat_pre(7);
    hxk_pre(14)=xhat_pre(18)-xhat_pre(9);    
    hxk_pre(15)=xhat_pre(19)-xhat_pre(11);
    
    
    jh_11=-g*cos(xhat_pre(1));
    
    jh_21=(-xhat_pre(6)*sin(xhat_pre(1))*sin(xhat_pre(3)));
    jh_22=cos(xhat_pre(3));                    
    jh_23=(-xhat_pre(2)*sin(xhat_pre(3)))...
         +(xhat_pre(6)*cos(xhat_pre(1))*cos(xhat_pre(3)));
    jh_26=cos(xhat_pre(1))*sin(xhat_pre(3));
    
    jh_31=-g*sin(xhat_pre(1))*sin(xhat_pre(3));
    jh_33=g*cos(xhat_pre(1))*cos(xhat_pre(3));
    
    jh_41=-xhat_pre(6)*cos(xhat_pre(1));
    jh_46=-sin(xhat_pre(1));     
    
    jh_51=(-xhat_pre(6)*sin(xhat_pre(1))*cos(xhat_pre(3)));
    jh_52=-sin(xhat_pre(3));
    jh_53=-(xhat_pre(2)*cos(xhat_pre(3)))...
         -(xhat_pre(6)*cos(xhat_pre(1))*sin(xhat_pre(3)));
    jh_56=cos(xhat_pre(1))*cos(xhat_pre(3));

    jh_61=(-xhat_pre(7)*sin(xhat_pre(1))*cos(xhat_pre(5)))...
         +(xhat_pre(9)*cos(xhat_pre(1))*sin(xhat_pre(3))*cos(xhat_pre(5)))...
         +(xhat_pre(11)*cos(xhat_pre(1))*cos(xhat_pre(3))*cos(xhat_pre(5)));
    jh_63=(xhat_pre(9)*(sin(xhat_pre(5))*sin(xhat_pre(3))...
                +sin(xhat_pre(1))*cos(xhat_pre(3))*cos(xhat_pre(5))) )...
         +(xhat_pre(11)*(cos(xhat_pre(3))*sin(xhat_pre(5))...
                -sin(xhat_pre(1))*sin(xhat_pre(3))*cos(xhat_pre(5))) );
    jh_65=(-xhat_pre(7)*cos(xhat_pre(1))*sin(xhat_pre(5)))...
         +(xhat_pre(9)*(-cos(xhat_pre(5))*cos(xhat_pre(3))...
                -sin(xhat_pre(1))*sin(xhat_pre(3))*sin(xhat_pre(5))) )...
         +(xhat_pre(11)*(sin(xhat_pre(3))*cos(xhat_pre(5))...
                -sin(xhat_pre(1))*cos(xhat_pre(3))*sin(xhat_pre(5))) );
    jh_67=cos(xhat_pre(1))*cos(xhat_pre(5));
    jh_69=-sin(xhat_pre(5))*cos(xhat_pre(3))...
                +sin(xhat_pre(1))*sin(xhat_pre(3))*cos(xhat_pre(5)) ;
    jh_611=sin(xhat_pre(3))*sin(xhat_pre(5))...
                +sin(xhat_pre(1))*cos(xhat_pre(3))*cos(xhat_pre(5));

    jh_71=(-xhat_pre(7)*sin(xhat_pre(5))*sin(xhat_pre(1)))...
         +(xhat_pre(9)*cos(xhat_pre(1))*sin(xhat_pre(3))*sin(xhat_pre(5)))...
         +(xhat_pre(11)*cos(xhat_pre(1))*cos(xhat_pre(3))*sin(xhat_pre(5)));
    jh_73=(xhat_pre(9)*(-sin(xhat_pre(3))*cos(xhat_pre(5))...
                +sin(xhat_pre(1))*cos(xhat_pre(3))*sin(xhat_pre(5))) )...
         +(xhat_pre(11)*(-cos(xhat_pre(3))*cos(xhat_pre(5))...
                -sin(xhat_pre(1))*sin(xhat_pre(3))*sin(xhat_pre(5))) );
    jh_75=(xhat_pre(7)*cos(xhat_pre(5))*cos(xhat_pre(1)))...
         +(xhat_pre(9)*(-cos(xhat_pre(3))*sin(xhat_pre(5))...
                +sin(xhat_pre(1))*sin(xhat_pre(3))*cos(xhat_pre(5))) )...
         +(xhat_pre(11)*(sin(xhat_pre(3))*sin(xhat_pre(5))...
                +sin(xhat_pre(1))*cos(xhat_pre(3))*cos(xhat_pre(5))) );
    jh_77=sin(xhat_pre(5))*cos(xhat_pre(1));
    jh_79=cos(xhat_pre(3))*cos(xhat_pre(5))...
                +sin(xhat_pre(1))*sin(xhat_pre(3))*sin(xhat_pre(5));
    jh_711=-sin(xhat_pre(3))*cos(xhat_pre(5))...
                +sin(xhat_pre(1))*cos(xhat_pre(3))*sin(xhat_pre(5));
    
    jh_81=(-xhat_pre(7)*cos(xhat_pre(1)))...
         -(xhat_pre(9)*sin(xhat_pre(1))*sin(xhat_pre(3)))...
         -(xhat_pre(11)*sin(xhat_pre(1))*cos(xhat_pre(3)));
    jh_83=(-xhat_pre(7)*sin(xhat_pre(1)))...
         +(xhat_pre(9)*cos(xhat_pre(1))*cos(xhat_pre(3)))...
         -(xhat_pre(11)*cos(xhat_pre(1))*sin(xhat_pre(3)));     
    jh_87=-sin(xhat_pre(1));
    jh_89=cos(xhat_pre(1))*sin(xhat_pre(3));
    jh_811=cos(xhat_pre(1))*cos(xhat_pre(3));
    
    jh_91=-B(1)*cos(xhat_pre(5))*sin(xhat_pre(1))...
         -B(2)*sin(xhat_pre(5))*sin(xhat_pre(1))...
         -B(3)*cos(xhat_pre(1));
    jh_95=-B(1)*sin(xhat_pre(5))*cos(xhat_pre(1))...
         +B(2)*cos(xhat_pre(5))*cos(xhat_pre(1));
     
    jh_101=B(1)*(cos(xhat_pre(5))*sin(xhat_pre(3))*cos(xhat_pre(1)))...
          +B(2)*(sin(xhat_pre(5))*sin(xhat_pre(3))*cos(xhat_pre(1)))...
          -B(3)*sin(xhat_pre(1))*sin(xhat_pre(3));    
    jh_103=B(1)*(sin(xhat_pre(5))*sin(xhat_pre(3))...
                    +cos(xhat_pre(5))*cos(xhat_pre(3))*sin(xhat_pre(1)))...
          +B(2)*(-cos(xhat_pre(5))*sin(xhat_pre(3))...
                    +sin(xhat_pre(5))*cos(xhat_pre(3))*sin(xhat_pre(1)))...
          +B(3)*cos(xhat_pre(1))*cos(xhat_pre(3));         
    jh_105=B(1)*(-cos(xhat_pre(5))*cos(xhat_pre(3))...
                    -sin(xhat_pre(5))*sin(xhat_pre(3))*sin(xhat_pre(1)))...
          +B(2)*(-sin(xhat_pre(5))*cos(xhat_pre(3))...
                    +cos(xhat_pre(5))*sin(xhat_pre(3))*sin(xhat_pre(1)));
           
    jh_111=B(1)*(cos(xhat_pre(5))*cos(xhat_pre(3))*cos(xhat_pre(1)))...
          +B(2)*(sin(xhat_pre(5))*cos(xhat_pre(3))*cos(xhat_pre(1)))...
          -B(3)*sin(xhat_pre(1))*cos(xhat_pre(3));               
    jh_113=B(1)*(sin(xhat_pre(5))*cos(xhat_pre(3))...
                    -cos(xhat_pre(5))*sin(xhat_pre(3))*sin(xhat_pre(1)))...
          +B(2)*(-cos(xhat_pre(5))*cos(xhat_pre(3))...
                    -sin(xhat_pre(5))*sin(xhat_pre(3))*sin(xhat_pre(1)))...
          -B(3)*cos(xhat_pre(1))*sin(xhat_pre(3));               
    jh_115=B(1)*(cos(xhat_pre(5))*sin(xhat_pre(3))...
                    -sin(xhat_pre(5))*cos(xhat_pre(3))*sin(xhat_pre(1)))...
          +B(2)*(sin(xhat_pre(5))*sin(xhat_pre(3))...
                    +cos(xhat_pre(5))*cos(xhat_pre(3))*sin(xhat_pre(1)));     
                
    jh_121=-g*sin(xhat_pre(1))*cos(xhat_pre(3));
    jh_123=-g*cos(xhat_pre(1))*sin(xhat_pre(3));

    J_H = ...
    [jh_11   0     0    0   0     0     0   1    0   0   0    0 0 0 0 0 0 0 0; %accX
     jh_21 jh_22 jh_23  0   0   jh_26   0   0    0   0   0    0 1 0 0 0 0 0 0; %gyroY
     jh_31   0   jh_33  0   0     0     0   0    0   1   0    0 0 0 0 0 0 0 0; %accY
     jh_41   0     0    1   0   jh_46   0   0    0   0   0    0 0 1 0 0 0 0 0; %gyroX
     jh_51 jh_52 jh_53  0   0   jh_56   0   0    0   0   0    0 0 0 1 0 0 0 0; %gyroZ
     jh_61   0   jh_63  0 jh_65   0   jh_67 0  jh_69 0 jh_611 0 0 0 0 0 0 0 0; %gpsx
     jh_71   0   jh_73  0 jh_75   0   jh_77 0  jh_79 0 jh_711 0 0 0 0 0 0 0 0; %gpsy
     jh_81   0   jh_83  0   0     0   jh_87 0  jh_89 0 jh_811 0 0 0 0 0 0 0 0; %gpsz
     jh_91   0     0    0 jh_95   0     0   0    0   0   0    0 0 0 0 0 0 0 0; %Bx
     jh_101  0   jh_103 0 jh_105  0     0   0    0   0   0    0 0 0 0 0 0 0 0; %By
     jh_111  0   jh_113 0 jh_115  0     0   0    0   0   0    0 0 0 0 0 0 0 0; %Bz 
     jh_121  0   jh_123 0   0     0     0   0    0   0   0    1 0 0 0 1 0 0 0; %accZ
       0     0     0    0   0     0    -1   0    0   0   0    0 0 0 0 0 1 0 0; %sonicu       
       0     0     0    0   0     0     0   0   -1   0   0    0 0 0 0 0 0 1 0; %sonicv
       0     0     0    0   0     0     0   0    0   0  -1    0 0 0 0 0 0 0 1];%sonicw

    S=(J_H*P_pre*J_H')+R;
    K=(P_pre*J_H')/S;                  
    innov=(z(:,i)-hxk_pre);
    update=K*innov;
    P_post=(eye(19)-K*J_H)*P_pre;            
    xhat(:,i)=xhat_pre+update;                 
end
end

