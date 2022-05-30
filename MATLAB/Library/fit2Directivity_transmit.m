function [n_output] = fit2Directivity_transmit(D_transmit)
%Calculates power of n based on directivity
% Inputs: D [dBi]
% Outputs: n

n=linspace(1,100,1001);

thetaB = sqrt(4*pi/(10^(D_transmit/10))); %desired HPBW

theta=linspace(0,pi/2,1001);

for nn=1:length(n)
    DR=cos(theta).^n(nn);
    HPBW(nn)=theta(find(DR<=sqrt(1/2),1,'first'));
end

n_output=n(find(HPBW<=thetaB,1));

end
