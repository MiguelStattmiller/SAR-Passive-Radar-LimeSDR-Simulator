function [Wi_ref] = RXGaussian_ref(X_receiver,Y_receiver,X_transmitter,Y_transmitter,alpha_zeroSVref,D_ref)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** INPUTS

alpha_zeroSVref=2 % Introduce desired receiver antenna angle in radians
D_ref=1 % Antenna length in meters


%**************** Function Calculus

% Define vector between receiver and transmitter
V=[X_receiver-X_transmitter Y_receiver-Y_transmitter]


u = [cos(alpha_zeroSVref),sin(alpha_zeroSVref)]; % direction beam of receiver antenna
Vectors_product=dot( V,u)/norm(V); % Calculate angle between vector target-receiver and direction beam antenna vector
DoA_ref =acos(Vectors_product);
tetha_3dBref= sqrt((4*pi)/D_ref);



Wi_ref=exp(-(DoA_ref-alpha_zeroSVref).^2 / (c*(tetha_3dBref.^2))) 