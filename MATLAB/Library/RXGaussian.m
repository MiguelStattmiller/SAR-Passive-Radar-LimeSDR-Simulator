function [Wi_Surv] = RXGaussian(X_receiver,Y_receiver,X_target,Y_target,alpha_zero_surv,D_surv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%**************** Function Calculus
theta_3dB_surv=sqrt(4*pi/(10^(D_surv/10))); %desired HPBW
n=fit2Directivity_surv(D_surv);
alpha=linspace(alpha_zero_surv-theta_3dB_surv/2,alpha_zero_surv+theta_3dB_surv/2);

% Define vector between receiver and target

V=[X_target-X_receiver Y_target-Y_receiver];
alpha_DoA =atan2(abs(V(1)), V(2));

%u = [cos(alpha_zero_surv),sin(alpha_zero_surv)]; % direction beam of receiver antenna
%alpha_zero = atan2(abs(u(2)), u(1));

if  and(alpha_DoA>alpha(1), alpha_DoA<alpha(end))
    Wi_Surv= cos(alpha_DoA-alpha_zero_surv)^n;
else
     Wi_Surv= 0;
end



end