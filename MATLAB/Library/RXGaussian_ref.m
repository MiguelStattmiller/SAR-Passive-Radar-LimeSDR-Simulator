function [Wi_ref] = RXGaussian_ref(X_receiverref,Y_receiverref,X_transmitter,Y_transmitter,alpha_zeroSVref,D_ref)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** INPUTS

alpha_zeroSVref=2 % Introduce desired receiver antenna angle in radians
D_ref=1 % Antenna length in meters
n=1; % Antenna directivity


%**************** Function Calculus

% Define vector between receiver and transmitter
V=[X_receiverref-X_transmitter Y_receiverref-Y_transmitter];
angle = atan2(abs(V(2)), V(1));


u = [cos(alpha_zeroSVref),sin(alpha_zeroSVref)]; % direction beam of receiver antenna
angle2 = atan2(abs(u(2)), u(1));




Wi_ref=cos(angle-angle2)^n;