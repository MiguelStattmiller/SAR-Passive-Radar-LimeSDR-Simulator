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
dev.rx1.frequency   = Fc;
dev.rx1.samplerate  = Fs;
dev.rx1.bandwidth   = BW;
dev.rx1.gain        = Gain;
dev.rx1.antenna     = 1;     % LNAH

% Read parameters from the devices
Fs_dev1      = dev.rx1.samplerate;  
Fc_dev1      = dev.rx1.frequency;
BW_dev1      = dev.rx1.bandwidth;
Ant_dev1     = dev.rx1.antenna;
Gain_dev1    = dev.rx1.gain;


fprintf('Rx Device sampling frequency: %3.1fHz, Initial sampling frequency: %3.1fHz\n', Fs_dev1, Fs);
fprintf('Rx Device carrier frequency: %3.1fHz, Initial carrier frequency: %3.1fHz\n', Fc_dev1, Fc);
fprintf('Rx Device bandwidth: %3.1fHz, Initial bandwith: %3.1fHz\n', BW_dev1, BW);
fprintf('Rx Device antenna: %d \n', Ant_dev1);
fprintf('Rx Device gain: %3.1fdB, Initial gain: %3.1fdB\n', Gain_dev1, Gain);

% Create empty array for the received signal

bufferRx1    = complex(zeros(TotalTime*Fs,1));

% Enable stream parameters. These may NOT be changed while the device is streaming.

dev.rx1.enable;


% Start the module
dev.start();
fprintf('Start of LimeSDR\n');
tic;
start_time=datetime('now','Format','dd-MMM-uuuu HH:mm:ss.SSS');

% Receive samples on RX1 channel
indRx1 = 1;  % index of the last received sample

    [samples1, ~, samplesLength1]             = dev.receive(Fs*Ts,1);
    bufferRx1(indRx1:indRx1+samplesLength1-1)   = samples1;

pause(1)


% Cleanup and shutdown by stopping the RX stream and having MATLAB delete the handle object.
dev.stop();
tempo_rececao=toc;
stop_time=datetime('now','Format','dd-MMM-uuuu HH:mm:ss.SSS');
clear dev;
fprintf('Stop of LimeSDR\n');


TT=array2timetable(samples1,'SampleRate',Fs,'StartTime',start_time);
writetimetable(TT,'Timetable_Rx1.txt','Delimiter','bar');
