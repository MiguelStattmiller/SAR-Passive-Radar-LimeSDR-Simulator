% Author of the LimeSDR MATLAB compatibility program:
%    Damir Rakhimov, CRL, TU Ilmenau, Dec 2019

% Author of the time2freq function:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022

%This program allows the detection of signals with SDR in two channels.
clc
clear all

addpath('C:\Program Files\PothosSDR') % add path with LimeSuite library

% Initialize parameters
TotalTime   = 1;       % Time of observation, s
Fc          = 2.45e9;   % Carrier Frequency, Hz
Fs          = 10e6;      % Frequency of sampling frequency, Hz   2x> B
Ts          = 1;      % Signal duration, s
Fsig        = 2.45e9;    % Frequency of desired signal, Hz
Asig        = 1;        % Amplitude of signal, V
BW          = 5e6;      % Bandwidth of the signal, Hz (5-40MHz and 50-130Mhz)
Gain        = 20;       % Receiver Gain, dB
tempo=0:1/Fs:10000000/Fs-1/Fs;
fontsize=12;

% Open LimeSDR:
dev = limeSDR(); % Open device

% Setup device parameters. These may be changed while the device is actively streaming.
dev.rx0.frequency   = Fc;
dev.rx0.samplerate  = Fs;
dev.rx0.bandwidth   = BW;
dev.rx0.gain        = Gain;
dev.rx0.antenna     = 1;     % LNAH

dev.rx1.frequency   = Fc;
dev.rx1.samplerate  = Fs;
dev.rx1.bandwidth   = BW;
dev.rx1.gain        = Gain;
dev.rx1.antenna     = 1;     % LNAH

% Read parameters from the devices
Fs_dev      = dev.rx0.samplerate;  
Fc_dev      = dev.rx0.frequency;
BW_dev      = dev.rx0.bandwidth;
Ant_dev     = dev.rx0.antenna;
Gain_dev    = dev.rx0.gain;

Fs_dev1      = dev.rx1.samplerate;  
Fc_dev1      = dev.rx1.frequency;
BW_dev1      = dev.rx1.bandwidth;
Ant_dev1     = dev.rx1.antenna;
Gain_dev1    = dev.rx1.gain;

fprintf('Rx Device sampling frequency: %3.1fHz, Initial sampling frequency: %3.1fHz\n', Fs_dev, Fs);
fprintf('Rx Device carrier frequency: %3.1fHz, Initial carrier frequency: %3.1fHz\n', Fc_dev, Fc);
fprintf('Rx Device bandwidth: %3.1fHz, Initial bandwith: %3.1fHz\n', BW_dev, BW);
fprintf('Rx Device antenna: %d \n', Ant_dev);
fprintf('Rx Device gain: %3.1fdB, Initial gain: %3.1fdB\n', Gain_dev, Gain);


% Create empty array for the received signal
bufferRx    = complex(zeros(TotalTime*Fs,1));
bufferRx1    = complex(zeros(TotalTime*Fs,1));

% Enable stream parameters. These may NOT be changed while the device is streaming.
dev.rx0.enable;
dev.rx1.enable;


% Start the module
dev.start();
fprintf('Start of LimeSDR\n');


% Receive samples on RX0 channel
indRx = 1;
    [samples, ~, samplesLength]             = dev.receive(Fs*Ts,0);
    bufferRx(indRx:indRx+samplesLength-1)   = samples;


% Receive samples on RX1 channel
indRx1 = 1;  % index of the last received sample

    [samples1, ~, samplesLength1]             = dev.receive(Fs*Ts,1);
    bufferRx1(indRx1:indRx1+samplesLength1-1)   = samples1;

pause(1)

% Cleanup and shutdown by stopping the RX stream and having MATLAB delete the handle object.
dev.stop();
clear dev;
fprintf('Stop of LimeSDR\n');


%**************** Plot spectrograms of the received signals

figure(1)
subplot(3,2,1);
spectrogram(bufferRx,2^12,2^10,2^12,'centered','yaxis')
title('Reference Signal');

subplot(3,2,2);
spectrogram(bufferRx1,2^12,2^10,2^12,'centered','yaxis')
title('Surveillance Signal');

%**************** Plot Received signals Time Domaim

fig=figure;
set(fig,'color','white');
plot(tempo,20*log10(bufferRx),'linewidth',2,'color','b');
xlabel('Time [s]');
ylabel('Reference signal [dB]');
set(gca,'fontsize',fontsize);
grid on;

fig=figure;
set(fig,'color','white');
plot(tempo,20*log10(bufferRx1),'linewidth',2,'color','b');
xlabel('Time [s]');
ylabel('Surveillance signal [dB]');
set(gca,'fontsize',fontsize);
grid on;


%**************** Plot Received signals Frequency Domain

[Rs_FreqD,Rs_FD]=time2freq(bufferRx,tempo);
fig=figure;
set(fig,'color','white');
plot(abs(Rs_FreqD),20*log10(abs(Rs_FD)),'linewidth',2,'color','b');
xlabel('Frequency [Hz]');
ylabel('Reference signal [dB]');
set(gca,'fontsize',fontsize);
grid on;


[Rs_FreqD,Rs_FD]=time2freq(bufferRx1,tempo);
fig=figure;
set(fig,'color','white');
plot(abs(Rs_FreqD),20*log10(abs(Rs_FD)),'linewidth',2,'color','b');
xlabel('Frequency [Hz]');
ylabel('Surveillance signal [dB]');
set(gca,'fontsize',fontsize);
grid on;


% Select a few samples to get the process quicker
t = bufferRx(1:1000);
x = transpose(t);


t1 = bufferRx1(1:1000);
x1 = transpose(t1);


%**************** Frequency Domain of the receiving Samples

[freq,Spectrum]=time2freq(x,tempo);
fig=figure;
hold on
set(fig,'color','white');
plot(abs(freq),20*log10(abs(Spectrum)),'b','linewidth',2);
xlabel('Frequency [Hz]');
xlim auto; 
ylabel('Spectrum');
set(gca,'fontsize',fontsize);
grid on;
hold off


[freq1,Spectrum1]=time2freq(x1,tempo);
fig=figure;
hold on
set(fig,'color','white');
plot(abs(freq1),20*log10(abs(Spectrum1)),'b','linewidth',2);
xlabel('Frequency [Hz]');
ylabel('Spectrum');
set(gca,'fontsize',fontsize);
grid on;
hold off

%**************** Plot the receiving Samples Time Domain

fig=figure;
set(fig,'color','white');
plot(tempo,20*log10(abs(x)),'linewidth',2,'color','b');
xlabel('Time [s]');
ylabel('Reference signal [dB]');
set(gca,'fontsize',fontsize);
grid on;

fig=figure;
set(fig,'color','white');
plot(tempo,20*log10(abs(x1)),'linewidth',2,'color','b');
xlabel('Time [s]');
ylabel('Surveillance signal [dB]');
set(gca,'fontsize',fontsize);
grid on;


%**************** Calculate ambiguity and cross-ambiguity functions

%Sref ambiguity function
[afmag,delay] = ambgfun(bufferRx,Fs,250000,'cut','Doppler');
afmag = afmag*1; % Select plot gain *1
afmag(afmag>1 )= 1;

 %Sr ambiguity function
[afmag2,delay2] = ambgfun(bufferRx1,Fs,250000,'cut','Doppler');
afmag2 = afmag2*1; % Select plot gain *1
afmag2(afmag2>1 )= 1;

% Sref and Sr Cross-ambiguity function
[afmag3,delay3] = ambgfun(bufferRx,bufferRx1,Fs,[250000 250000],'cut','Doppler');
afmag3 = afmag3*1; % Select plot gain *1
afmag3(afmag3>1 )= 1;



%Plot the ambiguity and cross-ambiguity functions of Sref and Sr

%Ambiguity Function Sref
[pks1, index] = max(afmag);
xMax = delay(index);
yMax = pks1;
subplot(3,2,1)
plot(delay,afmag,'LineStyle','-.');
hold on
plot3(delay(pks1,index),afmag(pks1,index), pks1, '^r')
textString = sprintf('(%f, %f)', xMax, yMax);
text(xMax, yMax,textString,"Color",'b','FontSize',10);
hold off
shading interp;
xlim([-0.5e-5 0.5e-5]);
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');



%Ambiguity Function Sr
[pks2,index2] = max(afmag2);
xMax2 = delay2(index2);
yMax2 = pks2;
subplot(3,2,2)   
plot(delay2,afmag2,'LineStyle','-'); 
hold on
plot3(delay2(pks2,index2), afmag2(pks2,index2), pks2, '^r')
textString2 = sprintf('(%f, %f)', xMax2, yMax2);
text(xMax2, yMax2,textString2,"Color",'b','FontSize',10);
hold off
shading interp;
xlim([-0.5e-5 0.5e-5]);
grid on; 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sr');


% Plot the cross-ambiguity function of Sref and Sr
[pks3,index3] = max(afmag3);
xMax3 = delay3(index3);
yMax3 = pks3;
subplot(3,2,3)
plot(delay3,afmag3,'LineStyle','-'); 
hold on
%plot3(delay3(pks3,index3), afmag3(pks3,index3), pks3, '^r');
textString3 = sprintf('(%f, %f)', xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
%xlim([-0.5e-5 16e-5]);
xlim auto;
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Cross-correlation');

% Plot of Sref, Sr
subplot(3,2,4)
plot(delay,afmag,'LineStyle','-.','Color','g'); % Green Sref
hold on
plot(delay2, afmag2,'LineStyle','-','Color','r'); % Red Sr
plot(delay3,afmag3,'LineStyle','-','Color','b'); 
hold off
xlim ([-0.2 0.2]);
%xlim auto
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Sref, Sr and cross-ambiguity');
legend('Sref','Sr','Cross-ambiguity');

%************* Calculate ambiguity and cross-ambiguity functions-time and doppler delay

%Reference_Signal ambiguity function
[afmag,delay,doppler] = ambgfun(x,Fs,250000);
afmag = afmag*1; % Select plot gain *1
afmag(afmag>1 )= 1;

 %Surveillance_Signal ambiguity function
[afmag2,delay2,doppler2] = ambgfun(x1,Fs,250000);
afmag2 = afmag2*1; % Select plot gain *1
afmag2(afmag2>1 )= 1;

%Cross-ambiguity
[afmag3,delay3,doppler3] = ambgfun(x,x1,Fs,[250000 250000]);
afmag3 = afmag3*1; % Select plot gain *1
afmag3(afmag3>1 )= 1;

%**************** Plot ambiguity and cross-ambiguity functions time and doppler delay

%Plot Ambiguity Function of Sref
subplot(3,2,1)
surf(delay,doppler,afmag,'LineStyle','-.');
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');


%Plot Ambiguity Function of Sr
subplot(3,2,2)
surf(delay2,doppler2,afmag2,'LineStyle','-.');
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');


% Plot cross-ambiguity function of Sref and Sr

figure;
subplot(3,2,3)
surf(delay3,doppler3,afmag3,'LineStyle','-.');
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');

