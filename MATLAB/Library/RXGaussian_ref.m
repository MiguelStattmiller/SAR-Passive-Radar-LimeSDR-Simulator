function [Wi_ref] = RXGaussian_ref(X_receiverref,Y_receiverref,X_transmitter,Y_transmitter,alpha_zeroSVref,theta_3dB_ref)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%**************** Function Calculus
n=real((1/2*log10(2))/(log10(cos(theta_3dB_ref))));
alpha=linspace(alpha_zeroSVref-pi/2,pi/2+alpha_zeroSVref);

% Define vector between receiver and transmitter
V=[X_receiverref-X_transmitter Y_receiverref-Y_transmitter];
alpha_DoA = atan2(abs(V(2)), V(1));


%u = [cos(alpha_zeroSVref),sin(alpha_zeroSVref)]; % direction beam of receiver antenna
%angle2 = atan2(abs(u(2)), u(1));

if  and(alpha_DoA>alpha(1), alpha_DoA<alpha(end))
    Wi_ref=cos(alpha_DoA-alpha_zeroSVref)^n;
else
     Wi_ref= 0;
end


