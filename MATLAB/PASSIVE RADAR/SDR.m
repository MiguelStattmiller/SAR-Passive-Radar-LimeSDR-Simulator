
% Author of the LimeSDR MATLAB compatibility program:
%    Damir Rakhimov, CRL, TU Ilmenau, Dec 2019

% Author of the current program, based on the simple RX by Damir Rakhimov:
% Afonso Sénica, Escola Naval, June 2020
clc
clear all

addpath('C:\Program Files\PothosSDR\bin') % add path with LimeSuite library

% Initialize parameters
TotalTime   = 12;       % Time of observation, s
Fc          = 1003e6;   % Carrier Frequency, Hz
Fs          = 1e6;      % Frequency of sampling frequency, Hz
Ts          = 4e0;      % Signal duration, s
Fsig        = 0.5e0;    % Frequency of desired signal, Hz
Fdev        = 1e5;      % Frequency of deviation, Hz
Fif         = 2e5;      % Intermediate frequency, Hz
Asig        = 1;        % Amplitude of signal, V
BW          = 5e6;      % Bandwidth of the signal, Hz (5-40MHz and 50-130Mhz)
Gain        = 10;       % Receiver Gain, dB

% Open LimeSDR:
dev = limeSDR(); % Open device

% Setup device parameters. These may be changed while the device is actively streaming.
dev.rx0.frequency   = Fc;
dev.rx0.samplerate  = Fs;
dev.rx0.bandwidth   = BW;
dev.rx0.gain        = Gain;
dev.rx0.antenna     = 2;     % LNA_L

dev.rx1.frequency   = Fc;
dev.rx1.samplerate  = Fs;
dev.rx1.bandwidth   = BW;
dev.rx1.gain        = Gain;
dev.rx1.antenna     = 2;     % LNA_L

% Read parameters from the devices
Fs_dev      = dev.rx0.samplerate;  % in SPS
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

% Select a few samples to get the process quicker
t = bufferRx(1:1000);
x = transpose(t);


t1 = bufferRx1(1:1000);
x1 = transpose(t1);

 %Select plot gain *1
%Sref ambiguity function
[afmag,delay,doppler] = ambgfun(x,Fs,250000);
afmag = afmag*1;
afmag(afmag>1 )= 1;

 %Sr ambiguity function
[afmag2,delay2,doppler2] = ambgfun(x1,Fs,250000);
afmag2 = afmag2*1;
afmag2(afmag2>1 )= 1;

%Correlation
[afmag3,delay3,doppler3] = ambgfun(x,x1,Fs,[250000 250000]);
afmag3 = afmag3*1.5;
afmag3(afmag3>1 )= 1;


% Plot spectrograms of the recieved signals
figure(1)
subplot(3,2,1);
spectrogram(bufferRx,2^12,2^10,2^12,'centered','yaxis')

subplot(3,2,2);
spectrogram(bufferRx1,2^12,2^10,2^12,'centered','yaxis')


% Plot the ambiguity functions of Sref and Sr

subplot(3,2,3)
surf(delay,doppler,afmag,'LineStyle','none'); 
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Doppler f_d (Hz)');
title('Ambiguity Function Sref');


subplot(3,2,4)   
surf(delay2,doppler2,afmag2,'LineStyle','none'); 
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Doppler f_d (kHz)');
title('Ambiguity Function Sr');


% Plot the correlation of Sref and Sr

subplot(3,2,5)
surf(delay3,doppler3,afmag3,'LineStyle','none'); 
shading interp;
axis([-0.3e-6 0.3e-6 -1500 1500]); 
zlim([0 1]);
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Doppler f_d (Hz)');
title('Cross-correlation');

fprintf('Time for visualisation: %g\n', toc);
