
% This code allows to represent the FFT of a rectangular Pulse.


% Define the wave
t=[-3:0.01:3]; 
x = rectpuls(t);

%Plot
subplot(2,1,1); 
plot(t,x); 
title("Rectangular Function"); 
xlabel('Time Domain'); 
%FFT
X=fft(x); 
%Plot
subplot(2, 1, 2);
plot(t, fftshift(abs(X))); 
xlabel('Frequency Domain'); 
