% Author of the LimeSDR MATLAB compatibility program:
%    Damir Rakhimov, CRL, TU Ilmenau, Dec 2019

% Author of the time2freq function:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the ambiguity function program, based on the simple RX by Damir Rakhimov:
% Afonso Sénica, Escola Naval, June 2020

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022

clc
clear all

addpath('C:\Program Files\PothosSDR\bin') % add path with LimeSuite library

% Initialize parameters
TotalTime   = 1;       % Time of observation, s
Fc          = 2.45e9;   % Carrier Frequency, Hz
Fs          = 17e6;      % Frequency of sampling frequency, Hz   2x> B
Ts          = 1;      % Signal duration, s
Fsig        = 2.45e9;    % Frequency of desired signal, Hz
Asig        = 1;        % Amplitude of signal, V
BW          = 8e6;      % Bandwidth of the signal, Hz (5-40MHz and 50-130Mhz)
Gain        = 20;       % Receiver Gain, dB
tempo=0:1/Fs:1000/Fs-1/Fs;
fontsize=12;

% Open LimeSDR:
dev = limeSDR(); % Open device

% Setup device parameters. These may be changed while the device is actively streaming.
dev.rx0.frequency   = Fc;
dev.rx0.samplerate  = Fs;
dev.rx0.bandwidth   = BW;
dev.rx0.gain        = Gain;
dev.rx0.antenna     = 2;     % LNAL

dev.rx1.frequency   = Fc;
dev.rx1.samplerate  = Fs;
dev.rx1.bandwidth   = BW;
dev.rx1.gain        = Gain;
dev.rx1.antenna     = 2;     % LNAL

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

% Select a few samples to get the process quicker
t = bufferRx(1:1000);
x = transpose(t);


t1 = bufferRx1(1:1000);
x1 = transpose(t1);


%% Receiving channels Spectrum

[freq,Spectrum]=time2freq(x,tempo);
fig=figure;
hold on
set(fig,'color','white');
plot(freq,20*log10(abs(Spectrum)),'b','linewidth',2);
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
plot(freq1,20*log10(abs(Spectrum1)),'b','linewidth',2);
xlabel('Frequency [Hz]');
ylabel('Spectrum');
set(gca,'fontsize',fontsize);
grid on;
hold off




%% Select plot gain *1
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
afmag3 = afmag3*1;
afmag3(afmag3>1 )= 1;

%%
% Plot spectrograms of the recieved signals
figure(1)
subplot(3,2,1);
spectrogram(bufferRx,2^12,2^10,2^12,'centered','yaxis')

subplot(3,2,2);
spectrogram(bufferRx1,2^12,2^10,2^12,'centered','yaxis')


%% Plot the ambiguity functions of Sref and Sr

[maxValue1] = max(afmag(:));
subplot(3,2,1)
surf(delay,doppler,afmag,'LineStyle','none');
text(-0.5e-5,-0.5e-5,maxValue1,['\leftarrow Máximo = ' num2str(maxValue1)],'color','b','FontSize',10);
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Doppler f_d (Hz)');
title('Ambiguity Function Sref');


%%
[maxValue2] = max(afmag2(:));
subplot(3,2,2)   
surf(delay2,doppler2,afmag2,'LineStyle','none');
text(-0.5e-5,-0.5e-5,maxValue2,['\leftarrow Máximo = ' num2str(maxValue2)],'color','b','FontSize',10);
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Doppler f_d (kHz)');
title('Ambiguity Function Sr');


%%
% Plot the correlation of Sref and Sr

[maxValue3] = max(afmag3(:));
subplot(3,2,3)
surf(delay3,doppler3,afmag3,'LineStyle','none'); 
text(-0.5e-5,-0.5e-5,maxValue3,['\leftarrow Máximo = ' num2str(maxValue3)],'Color','b','FontSize',10);
shading interp;
axis([-0.5e-5 0.5e-5 -10000 10000]); 
zlim([0 1]);
grid on; 
view([140,35]); 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Doppler f_d (Hz)');
title('Cross-correlation');
Idx = find((afmag3(:) == max(afmag3(:)))) %Descobrir índice do máximo
[af3maxRow,af3maxCol] = ind2sub(size(afmag3), Idx) %Descobrir índice do máximo