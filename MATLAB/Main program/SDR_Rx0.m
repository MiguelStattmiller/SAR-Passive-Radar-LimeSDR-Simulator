clc
clear all

addpath('C:\Users\black\OneDrive\Ambiente de Trabalho\Tese\programação\SDR') % add path with LimeSuite library

% Initialize parameters
filename='Canal câmara 3';
TotalTime   = 1;       % Time of observation, s
Fc          = 2.40e9;   % Carrier Frequency, Hz
Fs          = 10e6;      % Frequency of sampling frequency, Hz   2x> B
Ts          = 10;      % Signal duration, s
Fsig        = 2.40e9;    % Frequency of desired signal, Hz
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

% Read parameters from the devices
Fs_dev      = dev.rx0.samplerate;  
Fc_dev      = dev.rx0.frequency;
BW_dev      = dev.rx0.bandwidth;
Ant_dev     = dev.rx0.antenna;
Gain_dev    = dev.rx0.gain;

fprintf('Rx Device sampling frequency: %3.1fHz, Initial sampling frequency: %3.1fHz\n', Fs_dev, Fs);
fprintf('Rx Device carrier frequency: %3.1fHz, Initial carrier frequency: %3.1fHz\n', Fc_dev, Fc);
fprintf('Rx Device bandwidth: %3.1fHz, Initial bandwith: %3.1fHz\n', BW_dev, BW);
fprintf('Rx Device antenna: %d \n', Ant_dev);
fprintf('Rx Device gain: %3.1fdB, Initial gain: %3.1fdB\n', Gain_dev, Gain);


% Create empty array for the received signal
bufferRx    = complex(zeros(TotalTime*Fs,1));

% Enable stream parameters. These may NOT be changed while the device is streaming.
dev.rx0.enable;

% Start the module
dev.start();
fprintf('Start of LimeSDR\n');
start_time=datetime('now','Format','dd-MMM-uuuu HH:mm:ss.SSS');

    
counter=1;
samples1=[];

while true


timestamp(counter,:)=datestr(now,'HH-MM-SS_FFF');

% Receive samples on RX0 channel
indRx = 1;
    [samples1_temp, ~, samplesLength]             = dev.receive(Fs*Ts,0);
    %bufferRx(indRx:indRx+samplesLength-1)   = samples1_temp;

samples1=[samples1 samples1_temp];
counter=counter+1

end

save(sprintf('%s_samples',filename),'samples1','-v7.3');
save(sprintf('%s_timestamp',filename),'timestamp','-v7.3');

% Cleanup and shutdown by stopping the RX stream and having MATLAB delete the handle object.
dev.stop();
stop_time=datetime('now','Format','dd-MMM-uuuu HH:mm:ss.SSS');
clear dev;
fprintf('Stop of LimeSDR\n');


