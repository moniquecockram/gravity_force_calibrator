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
% TorPeDo_Xbar=13.249;
% TorPeDo_Ybar=13.128;
wheel_radius=0.2;
syms TorPeDo_leng positive
% TorPeDo_leng=sqrt(2*0.3^2);
gamma=(acos((TorPeDo_leng^2-2*0.3^2)/(-2*0.3^2))-pi)/-2;
L=2;
l=0.5*(L-TorPeDo_leng*tan(gamma));
alpha=atan((0.3*cos(gamma)-wheel_radius)/l);
delta=atan((l*tan(alpha))/(0.6*sin(gamma/2)+l));

%gravity constant
G=6.67430E-11;

%define rotation period (s) and frequency (Hz)
T=1;

%define measurement time (s)
% syms t
t=[2:1E-3:3];
% t=0;

ux=[1 0]';
uy=[0 1]';

% progressbar(0.25);

%% For mx1

%as derived,define distance to point from masses as function of time (m)

r1_mx1=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+wheel_radius*cos((2*pi*t/T)+acos((1/(wheel_radius))*((l/(2*cos(alpha)))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))));
r2_mx1=0.5*(l/cos(alpha)+sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2))+wheel_radius*cos((2*pi*(t+T/2)/T)+acos((1/(wheel_radius))*(l/(2*cos(alpha))-0.5*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)))); %half cycle out

Beta_close=acos((-(2*wheel_radius)^2+(l/cos(alpha))^2+l^2+(TorPeDo_leng-l*tan(alpha))^2)/(2*(l/cos(alpha))*sqrt(l^2+(TorPeDo_leng-l*tan(alpha))^2)));
beta1_close=alpha+Beta_close/2-(Beta_close/2)*cos(2*pi*t/T);
beta2_close=alpha+Beta_close/2-(Beta_close/2)*cos(2*pi*(t+T/2)/T);

u1_mx1=ux*sqrt(1./(1+(tan(beta1_close)).^2))-uy*sqrt((tan(beta1_close)).^2./(1+(tan(beta1_close)).^2));
u2_mx1=ux*sqrt(1./(1+(tan(beta2_close)).^2))-uy*sqrt((tan(beta2_close)).^2./(1+(tan(beta2_close)).^2));

%calculating forces from m1 & m2 on point
F1_mx1=(G*TorPeDo_mx1*m1_wheel./(r1_mx1.^2)).*u1_mx1;
F2_mx1=(G*TorPeDo_mx1*m2_wheel./(r2_mx1.^2)).*u2_mx1;

%force from TorPeDo's arm holding m_x1 to the structure
rho_mx1=[cos(pi-gamma) -sin(pi-gamma);sin(pi-gamma) cos(pi-gamma)]*uy.*ones(2,size(t,2));

%force from other masses on TorPeDo
F3_mx1=(G*TorPeDo_mx1*TorPeDo_mx2./(TorPeDo_leng.^2))*-uy.*ones(2,size(t,2));
F4_mx1=(G*TorPeDo_mx1*TorPeDo_my2./((2*(L-l)).^2))*-ux.*ones(2,size(t,2));
F5_mx1=(G*TorPeDo_mx1*TorPeDo_my1./(0.6^2))*rho_mx1;

%% For mx2

%as derived,define distance to point from masses as function of time (m)
r2_mx2=r1_mx1;
r1_mx2=r2_mx1;

u1_mx2=ux*sqrt(1./(1+(tan(beta1_close)).^2))+uy*sqrt((tan(beta1_close)).^2./(1+(tan(beta1_close)).^2));
u2_mx2=ux*sqrt(1./(1+(tan(beta2_close)).^2))+uy*sqrt((tan(beta2_close)).^2./(1+(tan(beta2_close)).^2));

%calculating forces from m1 & m2 on point
F1_mx2=(G*TorPeDo_mx2*m1_wheel./(r1_mx1.^2)).*u1_mx2;
F2_mx2=(G*TorPeDo_mx2*m2_wheel./(r2_mx1.^2)).*u2_mx2;

%force from TorPeDo's arm holding m_x1 to the structure
rho_mx2=[cos(gamma) -sin(gamma);sin(gamma) cos(gamma)]*uy.*ones(2,size(t,2));

%force from other mass on TorPeDo
F3_mx2=(G*TorPeDo_mx1*TorPeDo_mx2./(TorPeDo_leng.^2))*uy.*ones(2,size(t,2));
F4_mx2=(G*TorPeDo_mx2*TorPeDo_my1./((2*(L-l)).^2))*-ux.*ones(2,size(t,2));
F5_mx2=(G*TorPeDo_mx2*TorPeDo_my2./(0.6^2))*rho_mx2;

%% For my1

%as derived,define distance to point from masses as function of time (m)
r1_my1=0.5*(sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha))^2)+sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha)+2*wheel_radius)^2))+wheel_radius*cos(2*pi*(t+T/2)/T+acos((1/(2*wheel_radius))*(sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha))^2)-sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha)+2*wheel_radius)^2))));
r2_my1=0.5*(sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha))^2)+sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha)+2*wheel_radius)^2))+wheel_radius*cos(2*pi*t/T+acos((1/(2*wheel_radius))*(sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha))^2)-sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha)+2*wheel_radius)^2))));

Beta_far=acos((2*(0.6*sin(gamma/2)+l)^2+(l*tan(alpha))^2+(l*tan(alpha)+2*wheel_radius)^2-(2*wheel_radius)^2)/(2*sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha))^2)*sqrt((0.6*sin(gamma/2)+l)^2+(l*tan(alpha)+2*wheel_radius)^2)));
beta1_far=delta+Beta_far/2-(Beta_far/2)*cos(2*pi*(t+T/2)/T);
beta2_far=delta+Beta_far/2-(Beta_far/2)*cos(2*pi*t/T);

u1_my1=ux*sqrt(1./(1+(tan(beta1_far)).^2))+uy*sqrt((tan(beta1_far)).^2./(1+(tan(beta1_far)).^2));
u2_my1=ux*sqrt(1./(1+(tan(beta2_far)).^2))+uy*sqrt((tan(beta2_far)).^2./(1+(tan(beta2_far)).^2));

%calculating forces from m1 & m2 on point
F1_my1=(G*TorPeDo_my1*m1_wheel./(r1_my1.^2)).*u1_my1;
F2_my1=(G*TorPeDo_my1*m2_wheel./(r2_my1.^2)).*u2_my1;

%force from TorPeDo's arm holding m_x1 to the structure
rho_my1=[cos(2*pi-gamma) -sin(2*pi-gamma);sin(2*pi-gamma) cos(2*pi-gamma)]*uy.*ones(2,size(t,2));

%force from other mass on TorPeDo
F3_my1=(G*TorPeDo_my1*TorPeDo_my2./(TorPeDo_leng.^2))*uy.*ones(2,size(t,2));
F4_my1=(G*TorPeDo_mx2*TorPeDo_my1./((2*(L-l)).^2))*ux.*ones(2,size(t,2));
F5_my1=(G*TorPeDo_mx1*TorPeDo_my1./(0.6^2))*rho_my1;

%% For my2

%as derived,define distance to point from masses as function of time (m)
r2_my2=r1_my1;
r1_my2=r2_my1;

u1_my2=ux*sqrt(1./(1+(tan(beta2_far)).^2))-uy*sqrt((tan(beta2_far)).^2./(1+(tan(beta2_far)).^2));
u2_my2=ux*sqrt(1./(1+(tan(beta1_far)).^2))-uy*sqrt((tan(beta1_far)).^2./(1+(tan(beta1_far)).^2));

%calculating forces from m1 & m2 on point
F1_my2=(G*TorPeDo_my2*m1_wheel./(r1_my2.^2)).*u1_my2;
F2_my2=(G*TorPeDo_my2*m2_wheel./(r2_my2.^2)).*u2_my2;

%force from TorPeDo's arm holding m_x1 to the structure
rho_my2=[cos(pi+gamma) -sin(pi+gamma);sin(pi+gamma) cos(pi+gamma)]*uy.*ones(2,size(t,2));

%force from other mass on TorPeDo
F3_my2=(G*TorPeDo_my1*TorPeDo_my2./(TorPeDo_leng.^2))*-uy.*ones(2,size(t,2));
F4_my2=(G*TorPeDo_mx1*TorPeDo_my2./((2*(L-l)).^2))*ux.*ones(2,size(t,2));
F5_my2=(G*TorPeDo_my1*TorPeDo_my2./(0.6^2))*rho_my2;

%% For All TorPeDo Masses

for i=1:size(t,2)
    A_mx2(:,i)=(dot(F2_mx2(:,i),rho_mx2(:,i))/norm(rho_mx2(:,i))).*rho_mx2(:,i);
    B_mx2(:,i)=(dot(F1_mx2(:,i),rho_mx2(:,i))/norm(rho_mx2(:,i))).*rho_mx2(:,i);
    A_mx1(:,i)=(dot(F2_mx1(:,i),-rho_mx1(:,i))/norm(-rho_mx1(:,i))).*-rho_mx1(:,i);
    B_mx1(:,i)=(dot(F1_mx1(:,i),-rho_mx1(:,i))/norm(-rho_mx1(:,i))).*-rho_mx1(:,i);
    A_my2(:,i)=(dot(F2_my2(:,i),rho_my2(:,i))/norm(rho_my2(:,i))).*rho_my2(:,i);
    B_my2(:,i)=(dot(F1_my2(:,i),rho_my2(:,i))/norm(rho_my2(:,i))).*rho_my2(:,i);
    A_my1(:,i)=(dot(F2_my1(:,i),-rho_my1(:,i))/norm(-rho_my1(:,i))).*-rho_my1(:,i);
    B_my1(:,i)=(dot(F1_my1(:,i),-rho_my1(:,i))/norm(-rho_my1(:,i))).*-rho_my1(:,i);
end

C_mx1=(dot(F3_mx1,rho_mx1)/norm(rho_mx1)).*rho_mx1;
D_mx1=(dot(F4_mx1,rho_mx1)/norm(rho_mx1)).*rho_mx1;
F6_mx1=A_mx1+B_mx1+C_mx1+D_mx1;

%force magnitude from both wheel masses on m_x1
F_tot_mx1=F1_mx1+F2_mx1+F3_mx1+F4_mx1-F5_mx1-F6_mx1;

C_mx2=(dot(F3_mx2,rho_mx2)/norm(rho_mx2)).*rho_mx2;
D_mx2=(dot(F4_mx2,rho_mx2)/norm(rho_mx2)).*rho_mx2;
F6_mx2=A_mx2+B_mx2+C_mx2+D_mx2;

%force magnitude from both wheel masses on m_x1
F_tot_mx2=F1_mx2+F2_mx2+F3_mx2+F4_mx2-F5_mx2-F6_mx2;

C_my1=(dot(F3_my1,rho_my1)/norm(rho_my1)).*rho_my1;
D_my1=(dot(F4_my1,rho_my1)/norm(rho_my1)).*rho_my1;
F6_my1=A_my1+B_my1+C_my1+D_my1;

%force magnitude from both wheel masses on m_x1
F_tot_my1=F1_my1+F2_my1+F3_my1+F4_my1-F5_my1-F6_my1;

C_my2=(dot(F3_my2,rho_my2)/norm(rho_my2)).*rho_my2;
D_my2=(dot(F4_my2,rho_my2)/norm(rho_my2)).*rho_my2;
F6_my2=A_my2+B_my2+C_my2+D_my2;

%force magnitude from both wheel masses on m_x1
F_tot_my2=F1_my2+F2_my2+F3_my2+F4_my2-F5_my2-F6_my2;

%% Combining wheel effects on TorPeDo
a=0;
b=0.6;
i=1;
d_T=zeros(1,size(t,2));
while i<size(t,2)
%     for k=1:6
    if (i==1)||(i==2)
        d_T(i)=double(vpasolve(F_tot_mx2(1,i)==F_tot_mx1(1,i),TorPeDo_leng,[a b]));
    elseif (d_T(i-1)==a) || (d_T(i-1)==b) || (d_T(i-1)>d_T(i-2)+0.03) || (d_T(i-1)<d_T(i-2)-0.03)
        d_T(i)=double(vpasolve(F_tot_mx2(2,i)==-F_tot_mx1(2,i),TorPeDo_leng,[a b]));
        d_T(i-1)=0.5*(d_T(i-4)+d_T(i-2));
    else
        d_T(i)=double(vpasolve(F_tot_mx2(2,i)==-F_tot_mx1(2,i),TorPeDo_leng,[a b]));
    end
    i=i+1;
%     end
end

for i=1:size(t,2)
    figure(1)
    hold on
    vectarrow([0 0],subs(F_tot_mx1(:,i),TorPeDo_leng,d_T(i)));
    figure(2)
    hold on
    vectarrow([0 0],subs(F_tot_mx2(:,i),TorPeDo_leng,d_T(i)));
    figure(3)
    hold on
    vectarrow([0 0],subs(F_tot_my1(:,i),TorPeDo_leng,d_T(i)));
    figure(4)
    hold on
    vectarrow([0 0],subs(F_tot_my2(:,i),TorPeDo_leng,d_T(i)));
    
hold on
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
plot(t(1:size(t,2)-1),d_T(1:size(t,2)-1))
ylabel('Cavity Length (m)')
yyaxis right
plot(t(1:size(t,2)-1),k(1:size(t,2)-1))
ylabel('Force Magnitude on TorPeDo Mass (N)')
xlabel('Time (s)')
title('TorPeDo Cavity Length Change from The Wheel')

% for i=1:size(t,2)
%     k2(i)=double(subs(norm(F1_mx1(:,i)),TorPeDo_leng,d_T(i)));
%     k3(i)=double(subs(norm(F1_mx2(:,i)),TorPeDo_leng,d_T(i)));
%     k4(i)=double(subs(norm(F2_mx1(:,i)),TorPeDo_leng,d_T(i)));
%     k5(i)=double(subs(norm(F2_mx2(:,i)),TorPeDo_leng,d_T(i)));
% end
% 
% figure(6)
% hold on
% plot(t,k2)
% plot(t,k3)
% plot(t,k4)
% plot(t,k5)
% legend('k1','k2','k3','k4')

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


% for i=1:size(t,2)
%     vectarrow([0 0],subs(u2_mx1(:,i),TorPeDo_leng,d_T2(i)))
% hold on
% end