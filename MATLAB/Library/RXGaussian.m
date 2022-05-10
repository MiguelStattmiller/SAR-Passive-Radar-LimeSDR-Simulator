function [Wi_Surv] = RXGaussian(X_receiver,Y_receiver,X_target,Y_target,alpha_zeroSV,D_surv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** INPUTS

alpha_zeroSV=2 % Introduce desired receiver antenna angle in radians
D_surv=1 % Antenna length in meters
n=1; % Antenna directivity

%**************** Function Calculus

% Define vector between receiver and target
V=[X_receiver-X_target Y_receiver-Y_target];
angle = atan2(abs(V(2)), V(1));

u = [cos(alpha_zeroSV),sin(alpha_zeroSV)]; % direction beam of receiver antenna
angle2 = atan2(abs(u(2)), u(1));



Wi_Surv= cos(angle-angle2)^n;












end