clc
clear all
Fs=10e6;

% Import Data

Reference_signal = load('Canal gabinete_samples.mat');
Reference_signal=cell2mat(struct2cell(Reference_signal));

Reference_time=load('Canal câmara_timestamp.mat');



Surveillance_signal=load('Canal 5_samples.mat');
Surveillance_signal=cell2mat(struct2cell(Surveillance_signal));

Surveillance_time=load('Canal 5_timestamp.mat');

% Reference Signal Time Base

val=['19-05-21_475'
     '19-05-33_280'];
StartEnd=datetime(val,'InputFormat','HH-mm-ss_SSS')
RS_time=0:1/Fs:numel(Reference_signal)/Fs;
val = string(val);
val = replace(val, '-', ':');
val = replace(val, '_', '.');
d = duration(val);
d.Format = 'hh:mm:ss.SSS';
Rdelay=seconds(diff(d)); % Difference in seconds


% Surveillance Signal Time Base

val2=['17-10-58_086'
    '17-11-09_923'];
StartEnd2=datetime(val2,'InputFormat','HH-mm-ss_SSS')
SS_time=0:1/Fs:numel(Surveillance_signal)/Fs;
val2 = string(val2);
val2 = replace(val2, '-', ':');
val2 = replace(val2, '_', '.');
d2 = duration(val2);
d2.Format = 'hh:mm:ss.SSS';
Sdelay=seconds(diff(d2)); % Difference in seconds

%**************** Processing time delay code of the received signals

nRef = numel(Reference_signal)/Fs;
nSurv = numel(Surveillance_signal)/Fs;
pulse_size = 500000;

 [freq_RS,RS]=time2freq(Reference_signal((kk-1)*pulse_size+(1:pulse_size)),RS_time((kk-1)*pulse_size+(1:pulse_size)));
 fig=figure;
 hold on
 set(fig,'color','white');
 plot(abs(freq_RS),20*log10(abs(RS)),'b','linewidth',2);
 %xlim ([20e6 40e6]);
 xlabel('Frequency [Hz]');
 ylabel('QPSK Spectrum [dB]');
 %set(gca,'fontsize',fontsize);
 grid on;
 hold off

  [freq_SS,SS]=time2freq(Surveillance_signal((kk-1)*pulse_size+(1:pulse_size)),SS_time((kk-1)*pulse_size+(1:pulse_size)));
 fig=figure;
 hold on
 set(fig,'color','white');
 plot(abs(freq_SS),20*log10(abs(SS)),'b','linewidth',2);
 %xlim ([20e6 40e6]);
 xlabel('Frequency [Hz]');
 ylabel('QPSK Spectrum [dB]');
 %set(gca,'fontsize',fontsize);
 grid on;
 hold off


 diferenca=20*log10(abs(SS))-20*log10(abs(RS));
 freq=abs(freq_SS)-abs(freq_RS);
 
% Time delay processing


for k = 1:nRef
    Reference_signal_samples =Reference_signal((k-1)*pulse_size+(1:pulse_size));

    for kk =1:nSurv

   
    Surveillance_signal_samples=Surveillance_signal((kk-1)*pulse_size+(1:pulse_size));


    [afmag3,delay3] = ambgfun(Reference_signal_samples,Surveillance_signal_samples,Fs,[250000 250000],'cut','Doppler');
    afmag3(afmag3>1 )= 1;  
    end
    
end


% Plot the cross-ambiguity function of Sref and Sr
[pks3,index3] = max(afmag3);
xMax3 = delay3(index3);
yMax3 = pks3;
subplot(3,2,3)
plot(delay3,afmag3,'LineStyle','-'); 
hold on
textString3 = sprintf('(%f, %f)', xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Cross-correlation');

% Doppler delay processing


for k = 1:nRef
    Reference_signal_samples =Reference_signal((k-1)*pulse_size+(1:pulse_size));

    for kk =1:nSurv

   
    Surveillance_signal_samples=Surveillance_signal((kk-1)*pulse_size+(1:pulse_size));



    [afmag3,doppler3] = ambgfun(Reference_signal_samples,Surveillance_signal_samples,Fs,[250000 250000],'Cut','Delay');
    afmag3(afmag3>1 )= 1;
    end
    
end

% Plot the cross-ambiguity function of Sref and Sr

[pks3,index3] = max(afmag3);
xMax3 = doppler3(index3);
yMax3 = pks3;
subplot(3,2,3)
plot(doppler3,afmag3,'LineStyle','-'); 
hold on
textString3 = sprintf('(%.2e, %f)', xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
grid on; 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Cross-ambiguity');


% Time and Doppler delay processing


for k = 1:nRef
    Reference_signal_samples =Reference_signal((k-1)*pulse_size+(1:pulse_size));

    for kk =1:nSurv

   
    Surveillance_signal_samples=Surveillance_signal((kk-1)*pulse_size+(1:pulse_size));
    [afmag3,delay3,doppler3] = ambgfun(Reference_signal_samples,Surveillance_signal_samples,Fs,[250000 250000]);
    afmag3(afmag3>1 )= 1;

    end
    
end



