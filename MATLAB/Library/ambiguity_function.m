clc
close all
clear all

% Read and info IQ file
info = audioinfo('noticiasrtl.wav');
[a,fs] = audioread('noticiasrtl.wav');

info2 = audioinfo('poprtl.wav');
[a2,fs2] = audioread('poprtl.wav');

% From 2 columns to a 1 column complex data
y = a;
y(:,2) = y(:,2)*1i;
y(:,1) = y(:,1) + y(:,2);
y(:,2) = [];

y2 = a2;
y2(:,2) = y2(:,2)*1i;
y2(:,1) = y2(:,1) + y2(:,2);
y2(:,2) = [];

% Column to vector - transpose
t = transpose(y);
x = t(1:1000);

t2 = transpose(y2);
x2 = t2(1:1000);


% Maths
[afmag,delay,doppler] = ambgfun(x,fs,10000);
afmag = afmag*2;
afmag(afmag>1 )= 1;

[afmag2,delay2,doppler2] = ambgfun(x2,fs2,10000);
afmag2 = afmag2*2;
afmag2(afmag2>1 )= 1;

% Figure plot
figure
subplot(2,2,1)   
plot(abs(fftshift(fft(t))));


subplot(2,2,2)
surf(delay,doppler,afmag,'LineStyle','none'); 
shading interp;
axis tight; 
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Doppler f_d (kHz)');
title('Ambiguity Function');

subplot(2,2,3)   
plot(abs(fftshift(fft(t2))));

subplot(2,2,4)   
surf(delay2,doppler2,afmag2,'LineStyle','none'); 
shading interp;
axis tight; 
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Doppler f_d (kHz)');
title('Ambiguity Function');