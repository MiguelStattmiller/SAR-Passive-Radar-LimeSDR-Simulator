function [Wi_transmit] = RXGaussian_transmit(X_target,Y_target,X_transmitter,Y_transmitter,alpha_zero_transmit,D_transmit)
%This function allows the calculation of a raised-cosine equation


%**************** Function Calculus
theta_3dB_transmit=sqrt(4*pi/(10^(D_transmit/10))); %desired HPBW
n=fit2Directivity_transmit(D_transmit);
alpha=linspace(alpha_zero_transmit-theta_3dB_transmit/2,alpha_zero_transmit+theta_3dB_transmit/2);

% Define vector between target and transmitter

V=[X_target-X_transmitter Y_target-Y_transmitter];
alpha_DoA = abs(atan2((V(1)), V(2)));

%u = [cos(alpha_zero_surv),sin(alpha_zero_surv)]; % direction beam of receiver antenna
%alpha_zero = atan2(abs(u(2)), u(1));

if  and(alpha_DoA>alpha(1), alpha_DoA<alpha(end))
    Wi_transmit= cos(alpha_DoA-alpha_zero_transmit)^n;
else
     Wi_transmit= 0;
end



end
