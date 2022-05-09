function [Wi_Surv] = RXGaussian(X_receiver,Y_receiver,X_target,Y_target,alpha_zeroSV,D_surv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** INPUTS

alpha_zeroSV=2 % Introduce desired receiver antenna angle in radians
D_surv=1 % Antenna length in meters


%**************** Function Calculus

% Define vector between receiver and pixel
V=[X_receiver-X_target Y_receiver-Y_target];


u = [cos(alpha_zeroSV),sin(alpha_zeroSV)]; % direction beam of receiver antenna
Vectors_product=dot( V,u)/norm(V); % Calculate angle between vector target-receiver and direction beam antenna vector
DoA =acos(Vectors_product);
tetha_3dB= sqrt((4*pi)/D_surv);



Wi_Surv=exp(-(DoA-alpha_zeroSV).^2 / (c*(tetha_3dB.^2))) 











end