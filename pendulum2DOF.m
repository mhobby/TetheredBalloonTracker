function [ xp ] = pendulum2DOF( ~, x )
%simulation to model anglualr displacement from vertical pendulum with 2DOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M J HOBBY (2013) mhobby1979@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
L=x(5); %tether length
a=x(6); %c/m where c is proportional to the damping visocosity of the air 
lift=x(7); %balloon lift
m=x(8); %mass of sonde + balloon + tether
kU=x(9); %balloon drag coefficient
U=x(10); %mean wind in U over sim period
g=9.81; %gravitational constant

Qth=(a*x(2));%/(cos(x(3))^2);
Qph=a*x(4);

torque_lift(1) = -(lift*sin(x(1)))/(m*L*cos(x(3))) ...
                    + (sign(U)*(kU*(U^2)*cos(x(1)))/(m*L*cos(x(3))));
torque_lift(2) = -(lift*sin(x(3)))/(m*L) ;

xp=zeros(10,1);
xp(1)=x(2);
xp(3)=x(4);
xp(2)=(( ((g/L)*sin(x(1)))+(2*x(2)*x(4)*sin(x(3))) ) / cos(x(3))) +torque_lift(1)- Qth;
xp(4)=((g/L)*cos(x(1))*sin(x(3)))-((x(2)^2)*sin(x(3))*cos(x(3))) +torque_lift(2)- Qph;
end