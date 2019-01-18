% createPendulumData.m
% create double pendulum simulation data
% Keisuke Fujii

clear ; close all
global fixed 

% initial parameters for pendulum 
phi1                = 1/8*pi;
dtphi1              = 0;
phi2                = 1/8*pi;
dtphi2              = 0;
g                   = 9.81; 
m1                  = 1; 
m2                  = 1; 
l1                  = 1; 
l2                  = 1;
duration            = 500; 
fps                 = 40;  
movie               = false;%true;

figure;

interval=[0, duration];
ivp=[phi1; dtphi1; phi2; dtphi2] ;
fixed = [g; m1; m2; l1; l2];

y = double_pendulum(ivp, duration, fps, movie);

Ang = y([1 3],1:end) ; 
save('./doublePendulum','Ang');

figure(2)
time = 1/fps:1/fps:duration ;
T = 800 ;
plot(time(1:T),Ang(1:2,1:T)')
legend('\theta_1','\theta_2')
xlabel('time (s)') ;
ylabel('angle (rad)');
box off
