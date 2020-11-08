clear
close all
clc

% progressbar(0);

%% Formalities
%assume wheel diameter is less than TorPeDo side length (see diagram)
%define constants (all SI units)
m1_wheel=5;
m2_wheel=5;
TorPeDo_mx1=5;
TorPeDo_mx2=5;
TorPeDo_my1=5;
TorPeDo_my2=5;
TorPeDo_Xbar=13.249;
TorPeDo_Ybar=13.128;
wheel_radius=0.2;
syms TorPeDo_leng positive
% TorPeDo_leng=sqrt(2*0.3^2);
gamma=(acos((TorPeDo_leng^2-2*0.3^2)/(-2*0.3^2))-pi)/-2;
L=2;
l=0.5*(L-TorPeDo_leng*tan(gamma));
alpha=atan((0.3*cos(gamma)-wheel_radius)/l);

%gravity constant
G=6.67430E-11;

%define rotation period (s) and frequency (Hz)
T=1;

%define measurement time (s)
% syms t
t=[1:1E-3:3];
% t=0;

ux=[1 0]';
uy=[0 1]';

% progressbar(0.25);

%% For mx1

%as derived,define distance to point from masses as function of time (m)
r1_mx1=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+2*wheel_radius*cos((2*pi*t/T)+acos(1/(2*wheel_radius)*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))));
r2_mx1=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+2*wheel_radius*cos((2*pi*(t+T/2)/T)+acos(1/(2*wheel_radius)*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)))); %half cycle out

Beta1_mx1=acos((-(2*wheel_radius)^2+(l/cos(alpha))^2+l^2+(TorPeDo_leng-l*tan(alpha))^2)/(2*(l/cos(alpha))*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)));
beta1_mx1=alpha+Beta1_mx1*sin(pi*t/T);
beta2_mx1=alpha+Beta1_mx1*sin(pi*(t+T/2)/T);

u1_mx1=ux*sqrt(1./(1+(tan(beta1_mx1)).^2))-uy*sqrt((tan(beta1_mx1)).^2./(1+(tan(beta1_mx1)).^2));
u2_mx1=ux*sqrt(1./(1+(tan(beta2_mx1)).^2))-uy*sqrt((tan(beta2_mx1)).^2./(1+(tan(beta2_mx1)).^2));

%calculating forces from m1 & m2 on point
F1_mx1=(G*TorPeDo_mx1*m1_wheel./(r1_mx1.^2)).*u1_mx1;
F2_mx1=(G*TorPeDo_mx1*m2_wheel./(r2_mx1.^2)).*u2_mx1;

%force from other mass on TorPeDo
F3_mx1=(G*TorPeDo_mx1*TorPeDo_mx2./(TorPeDo_leng.^2))*-uy.*ones(2,size(t,2));

%force from TorPeDo's arm holding m_x1 to the structure
rho_mx1=[cos(pi-gamma) -sin(pi-gamma);sin(pi-gamma) cos(pi-gamma)]*uy.*ones(2,size(t,2));

%% For mx2

%as derived,define distance to point from masses as function of time (m)
r2_mx2=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+2*wheel_radius*cos((2*pi*t/T)+acos(1/(2*wheel_radius)*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))));
r1_mx2=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+2*wheel_radius*cos((2*pi*(t+T/2)/T)+acos(1/(2*wheel_radius)*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)))); %half cycle out

Beta1_mx2=acos((-(2*wheel_radius)^2+(l/cos(alpha))^2+l^2+(TorPeDo_leng-l*tan(alpha))^2)/(2*(l/cos(alpha))*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)));
beta1_mx2=alpha+Beta1_mx2*sin(pi*t/T);
beta2_mx2=alpha+Beta1_mx2*sin(pi*(t+T/2)/T);

u1_mx2=ux*sqrt(1./(1+(tan(beta1_mx2)).^2))+uy*sqrt((tan(beta1_mx2)).^2./(1+(tan(beta1_mx2)).^2));
u2_mx2=ux*sqrt(1./(1+(tan(beta2_mx2)).^2))+uy*sqrt((tan(beta2_mx2)).^2./(1+(tan(beta2_mx2)).^2));

%calculating forces from m1 & m2 on point
F1_mx2=(G*TorPeDo_mx2*m1_wheel./(r1_mx1.^2)).*u1_mx2;
F2_mx2=(G*TorPeDo_mx2*m2_wheel./(r2_mx1.^2)).*u2_mx2;

%force from other mass on TorPeDo
F3_mx2=(G*TorPeDo_mx1*TorPeDo_mx2./(TorPeDo_leng.^2))*uy.*ones(2,size(t,2));

%force from TorPeDo's arm holding m_x1 to the structure
rho_mx2=[cos(gamma) -sin(gamma);sin(gamma) cos(gamma)]*uy.*ones(2,size(t,2));

%% For mx1 & mx2

for i=1:size(t,2)
    A_mx2(:,i)=(dot(F2_mx2(:,i),rho_mx2(:,i))/norm(rho_mx2(:,i))).*rho_mx2(:,i);
    B_mx2(:,i)=(dot(F1_mx2(:,i),rho_mx2(:,i))/norm(rho_mx2(:,i))).*rho_mx2(:,i);
    A_mx1(:,i)=(dot(F2_mx1(:,i),-rho_mx1(:,i))/norm(-rho_mx1(:,i))).*-rho_mx1(:,i);
    B_mx1(:,i)=(dot(F1_mx1(:,i),-rho_mx1(:,i))/norm(-rho_mx1(:,i))).*-rho_mx1(:,i);
end

C_mx1=(dot(F3_mx1,rho_mx1)/norm(rho_mx1)).*rho_mx1;
F4_mx1=A_mx1+B_mx1+C_mx1;

%force magnitude from both wheel masses on m_x1
F_tot_mx1=F1_mx1+F2_mx1+F3_mx1-F4_mx1;

C_mx2=(dot(F3_mx2,rho_mx2)/norm(rho_mx2)).*rho_mx2;
F4_mx2=A_mx2+B_mx2+C_mx2;

%force magnitude from both wheel masses on m_x1
F_tot_mx2=F1_mx2+F2_mx2+F3_mx2-F4_mx2;

%% Combining wheel effects on TorPeDo
a = -0.1;
b=0.5;
i=1;
d_T=zeros(1,size(t,2));
while i<size(t,2)
    if (i==1)
        d_T(i)=double(vpasolve(norm(F_tot_mx2(2,i))==norm(F_tot_mx1(2,i)),TorPeDo_leng,[a b]));
    elseif (d_T(i-1)==a) || (d_T(i-1)==b)
        d_T(i)=double(vpasolve(norm(F_tot_mx2(2,i))==norm(F_tot_mx1(2,i)),TorPeDo_leng,[a b]));
        d_T(i-1)=0.5*(d_T(i-2)+d_T(i));
    elseif (d_T(i-1)<0.3) || (d_T(i-1)>0.4)
        d_T(i)=0.3;
    else
        d_T(i)=double(vpasolve(norm(F_tot_mx2(2,i))==norm(F_tot_mx1(2,i)),TorPeDo_leng,[a b]));
    end
    i=i+1;
end

% plot(t(1:size(t,2)-1),d_T(1:size(t,2)-1))
% ylabel('Cavity Length')
% xlabel('t(ds)')
% title('TorPeDo Cavity Length Change from The Wheel')

for i=1:size(t,2)
    k(i)=double(subs(norm(F_tot_mx1(:,i)),TorPeDo_leng,d_T(i)));
end
figure(1)
hold on
yyaxis left
plot(t,d_T2)
ylabel('Cavity Length (m)')
yyaxis right
plot(t,k)
ylabel('Force Magnitude on TorPeDo Mass (N)')
xlabel('Time (s)')
title('TorPeDo Cavity Length Change from The Wheel')

% zlabel('|F_{tot}|')
% zlim([0 10])

%plot gravity fluctuations over time
% figure(1)
% fsurf(norm(F_tot_mx1))
% xlim([sqrt(2*0.3^2)-0.1 sqrt(2*0.3^2)+0.1])
% ylim([0 2]) %showing two cycles of wheel
% zlim([0 10E-2])
% % legend('closest','next closest','etc')
% xlabel('|d_T|')
% ylabel('t(ds)')
% zlabel('|F|')
% title('The Wheel''s Gravity Influence on m_{x1} Over Time')
%


%finding d_T as a function of F and t along d_T axis
% figure(2)
% options = optimoptions('solvername','Display','iter');
% d_T=solve(F_tot_mx1(2,1),TorPeDo_leng);
% fsurf(d_T)
