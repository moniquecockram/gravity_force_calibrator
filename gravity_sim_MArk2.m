clear
close all
clc

%% Formalities
%assume wheel diameter is less than TorPeDo side length (see diagram)
%define constants (all SI units)
m1_wheel=1000;
m2_wheel=1000;
TorPeDo_mx1=5;
TorPeDo_mx2=5;
TorPeDo_my1=5;
TorPeDo_my2=5;
TorPeDo_Xbar=13.249;
TorPeDo_Ybar=13.128;
wheel_radius=0.2;
syms TorPeDo_leng
gamma=(acos((TorPeDo_leng^2-2*0.3^2)/(-2*0.3^2))-pi)/-2;
L=1;
l=0.5*(L-TorPeDo_leng*tan(gamma));
alpha=atan((0.3*cos(gamma)-wheel_radius)/l);

%gravity constant
G=6.67430E-11;

%define rotation period (s) and frequency (Hz)
T=1;

%define measurement time (s)
% syms t
t=[0:0.1:2];

ux=[1 0]';
uy=[0 1]';

%% For m_x1

%as derived,define distance to point from masses as function of time (m)
r1=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+2*wheel_radius*cos((2*pi*t/T)+acos(1/(2*wheel_radius)*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))));
r2=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+2*wheel_radius*cos((2*pi*(t+T/2)/T)+acos(1/(2*wheel_radius)*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)))); %half cycle out

Beta1=-acos((-(2*wheel_radius)^2+(l/cos(alpha))^2+l^2+(TorPeDo_leng-l*tan(alpha))^2)/(2*(l/cos(alpha))*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)));
beta1=-alpha+Beta1*sin(pi*t/T);
Beta2=acos((-(2*wheel_radius)^2+(l/cos(alpha))^2+l^2+(TorPeDo_leng-l*tan(alpha))^2)/(2*(l/cos(alpha))*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)));
beta2=alpha+Beta2*sin(pi*(t+T/2)/T);

u1=ux*sqrt(1./(1+(tan(beta1)).^2))-uy*sqrt((tan(beta1)).^2./(1+(tan(beta1)).^2));
u2=ux*sqrt(1./(1+(tan(beta2)).^2))-uy*sqrt((tan(beta2)).^2./(1+(tan(beta2)).^2));

%calculating forces from m1 & m2 on point
F1=(G*TorPeDo_mx1*m1_wheel./(r1.^2)).*u1;
F2=(G*TorPeDo_mx1*m2_wheel./(r2.^2)).*u2;

%force from other mass on TorPeDo
F3=(G*TorPeDo_mx1*TorPeDo_mx2./(TorPeDo_leng.^2))*-uy;

%force from TorPeDo's arm holding m_x1 to the structure
rho=[cos(pi-gamma) -sin(pi-gamma);sin(pi-gamma) cos(pi-gamma)]*uy;

for i=1:size(t,2)
    A(:,i)=(dot(F2(:,i),-rho)/norm(-rho)).*-rho;
    B(:,i)=(dot(F1(:,i),-rho)/norm(-rho)).*-rho;
end

C=(dot(F3,rho)/norm(rho)).*rho;
F4=A+B+C;

%force magnitude from both wheel masses on m_x1
F_tot=F1+F2+F3-F4;

d_T=[sqrt(2*0.3^2)-0.01:0.2/10:sqrt(2*0.3^2)+0.01];
for i=1:size(t,2)
    for k=1:size(d_T,2)
        F(i,k)=double(subs(norm(F_tot(:,i)),TorPeDo_leng,d_T(k)));
    end
end

surf(d_T,t,F)
xlabel('|d_T|')
ylabel('t(ds)')
zlabel('|F_{tot}|')
zlim([0 1E-1])

%plot gravity fluctuations over time
figure(1)
fsurf(norm(F_tot))
xlim([sqrt(2*0.3^2)-0.1 sqrt(2*0.3^2)+0.1])
ylim([0 2]) %showing two cycles of wheel
zlim([0 10E-2])
% legend('closest','next closest','etc')
xlabel('|d_T|')
ylabel('t(ds)')
zlabel('|F|')
title('The Wheel''s Gravity Influence on m_{x1} Over Time')

%finding d_T as a function of F and t along d_T axis
figure(2)
options = optimoptions('solvername','Display','iter');
d_T=solve(F_tot(2,1),TorPeDo_leng);
fsurf(d_T)
