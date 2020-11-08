clear
close all
clc

%define masses (kg)
m1=1000;
m2=1000;

%gravity constant
G=6.67430E-11;

%define rotation period (s)
T=1;

%define rotation frequency (Hz)
f=1/T;

%define measurement time (s)
t=[0:0.01:1];

%define distance from COM of both masses (in starting config) (m)
d=[0.51:0.01:1];

hold on

for i=1:size(t,2)
    for k=1:size(d,2)
    %define distance to point from masses as function of time (m)
    r1=d(k)+0.5*cos(2*pi*f*t(i));
    r2=d(k)-0.5*cos(2*pi*f*t(i));
    
    %calculating forces from m1 & m2 on point
    F1(k)=G*m1*m2/(r1^2);
    F2(k)=G*m1*m2/(r2^2);
       
    %total force
    F(i,k)=F1(k)+F2(k);
    end   
end

%plot gravity fluctuations over time
figure(1)
surf(F)

% legend('closest','next closest','etc')
xlabel('d(cm)')
ylabel('t(ds)')
zlabel('|F|\times 10^{-11}(N)')
title('The Wheel''s Gravity Along d Over Time')

for i=1:size(t,2)
    for k=1:size(d,2)
    %define distance to point from masses as function of time (m)
    r1=d(k)+0.2*cos(2*pi*f*t(i));
    r2=d(k)-0.2*cos(2*pi*f*t(i));
    
    %calculating forces from m1 & m2 on point
    F1(k)=G*m1*m2/(r1^2);
    F2(k)=G*m1*m2/(r2^2);
       
    %total force
    F(i,k)=F1(k)+F2(k);
    end   
end

%plot gravity fluctuations over time
figure(2)
surf(F)

% legend('closest','next closest','etc')
xlabel('d(cm)')
ylabel('t(ds)')
zlabel('|F|\times 10^{-11}(N)')
title('The Wheel''s Gravity Along d Over Time')

for i=1:size(t,2)
    for k=1:size(d,2)
    %define distance to point from masses as function of time (m)
    r1=d(k)+0.7*cos(2*pi*f*t(i));
    r2=d(k)-0.7*cos(2*pi*f*t(i));
    
    %calculating forces from m1 & m2 on point
    F1(k)=G*m1*m2/(r1^2);
    F2(k)=G*m1*m2/(r2^2);
       
    %total force
    F(i,k)=F1(k)+F2(k);
    end   
end

%plot gravity fluctuations over time
figure(3)
surf(F)

% legend('closest','next closest','etc')
xlabel('d(cm)')
ylabel('t(ds)')
zlabel('|F|\times 10^{-11}(N)')
title('The Wheel''s Gravity Along d Over Time')
