function [Wi_Surv] = RXGaussian(X_receiver,Y_receiver,X_target,Y_target,alpha_zeroSV,D_surv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** INPUTS

alpha_zero_surv=3*pi/4; % Introduce desired receiver antenna angle in radians
D_surv=1; % Antenna length in meters
n=3;


%**************** Function Calculus
theta_3dB=sqrt(4*pi/D_surv);
alpha=linspace(alpha_zero_surv-pi/2,pi/2+alpha_zero_surv);

% Define vector between receiver and target

V=[X_receiver-X_target Y_receiver-Y_target];
alpha_DoA = atan2(abs(V(2)), V(1));

%u = [cos(alpha_zero_surv),sin(alpha_zero_surv)]; % direction beam of receiver antenna
%alpha_zero = atan2(abs(u(2)), u(1));

if  and(alpha_DoA>alpha(1), alpha_DoA<alpha(end))
    Wi_Surv= cos(alpha_DoA-alpha_zero_surv)^n;
else
     Wi_Surv= 0;
end



end