function [Wi_transmit] = RXGaussian_transmit(X_target,Y_target,X_transmitter,Y_transmitter,alpha_zero_transmit,D_transmit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** Function Calculus
theta_3dB_transmit=sqrt(4*pi/(10^(D_transmit/10))); %desired HPBW
n=fit2Directivity_transmit(D_transmit);
alpha=linspace(alpha_zero_transmit-pi/2,pi/2+alpha_zero_transmit);

% Define vector between target and transmitter

V=[X_target-X_transmitter Y_target-Y_transmitter];
alpha_DoA = atan2(abs(V(2)), V(1));

%u = [cos(alpha_zero_surv),sin(alpha_zero_surv)]; % direction beam of receiver antenna
%alpha_zero = atan2(abs(u(2)), u(1));

if  and(alpha_DoA>alpha(1), alpha_DoA<alpha(end))
    Wi_transmit= cos(alpha_DoA-alpha_zero_transmit)^n;
else
     Wi_transmit= 0;
end



end