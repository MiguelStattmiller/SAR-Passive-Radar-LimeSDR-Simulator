function [Wi_ref] = RXGaussian_ref(X_receiverref,Y_receiverref,X_transmitter,Y_transmitter,alpha_zeroSVref,D_ref)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%**************** Function Calculus
theta_3dB_ref=sqrt(4*pi/(10^(D_ref/10))); %desired HPBW
n=fit2Directivity_ref(D_ref);
alpha=linspace(alpha_zeroSVref-theta_3dB_ref/2,alpha_zeroSVref+theta_3dB_ref/2);

% Define vector between receiver and transmitter
V=[X_receiverref-X_transmitter Y_receiverref-Y_transmitter];
alpha_DoA =pi/2+atan2(abs(V(2)), V(1));


%u = [cos(alpha_zeroSVref),sin(alpha_zeroSVref)]; % direction beam of receiver antenna
%angle2 = atan2(abs(u(2)), u(1));

if  and(alpha_DoA>alpha(1), alpha_DoA<alpha(end))
    Wi_ref=cos(alpha_DoA-alpha_zeroSVref)^n;
else
     Wi_ref= 0;
end


