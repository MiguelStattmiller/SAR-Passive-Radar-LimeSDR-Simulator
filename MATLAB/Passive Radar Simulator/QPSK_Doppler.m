% Author of the time2freq function:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the QPSK modulator:
%     Joao Felicio

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022


% The current program is a simulator for testing doppler shifts with white noise addition, for moving targets
% using QPSK transmitted signal.
% QPSK signal is created based on a message written in .txt file
% For the passive radar: Reference signal= QPSK signal
%            tir            Surveillance signal=QPSK_signal delayed in Time 



%**************** QPSK Modulator 
clear
close all
clc

%***************** INPUTS ******************
fc=30e6; %Carrier frequency
c=3e8; % Speed of light 
Rb=fc/100; %Bitrate
%SNR=10; %SNR [dB]
fontsize=12;
SentenceFilename = 'Mensagem_short.txt';
MeasuredFilename = 'ImplantableAntenna2Reader.s2p';
%*******************************************

% Loads table with coding
[num,txt,table] = xlsread('TabelaCodificacao.xlsx') ;

% text2Transmit=
text2Transmit=read_txt(SentenceFilename);

%Converts characters to bits using coding table
sequence_str=[];
for cc=1:length(text2Transmit)
    char_str=text2Transmit(cc);
    for row=1:size(table,1)
        if strcmp(char_str,table(row,1))
            sequence_str=[sequence_str,char(table(row,3))];
            break;
        end
    end
end
%Converts to double
for cc=1:length(sequence_str)
    sequence(cc)=str2num(sequence_str(cc));
end

% Reference_signal uses this sequence of bits
sequence=[];
sequence=[1,0,0,1,0,0,1,0,0,1,1,1];
Nb=length(sequence); %Number of bits
Ns=Nb/2; %Number of symbols

%Carrier
Ac=1; %Amplitude
wc=2*pi*fc;
Tc=1/fc; %Period
Tb=1/Rb; %Bit period
fs=2*fc; %Sampling frequency
t_bit=linspace(0,Tb,round(Tb*fs)); %Time base of bit

Reference_Signal=[];
Binary_signal=[];
bit=[];
Eb=Ac^2*Tb; %Bit energy
S=Eb*Rb; %Power of QPSK signal

PolarRZ=2*sequence-1;
for ss=1:length(sequence)/2
    I=PolarRZ(2*ss-1)*Ac.*cos(wc*t_bit); %I-channel
    Q=PolarRZ(2*ss)*Ac.*sin(wc*t_bit); %Q-channel
    
    %Adding noise
    QPSK_temp=I+Q;
%     QPSK_temp = awgn(QPSK_temp,SNR,S);
    
    %Modulated QPSK signal
    Reference_Signal=[Reference_Signal,QPSK_temp];
end

t=linspace(0,Ns*Tb,length(Reference_Signal));


%**************** Plot QPSK signal over time
fig=figure;
set(fig,'color','white');
plot(t,Reference_Signal,'linewidth',2,'color','b');
xlabel('Time [s]');
ylabel('QPSK signal');
set(gca,'fontsize',fontsize);
%set(gca,'fontsize',fs);
grid on;

%**************** Plot frequency spectrum of QPSK signal

 [freq,Spectrum]=time2freq(Reference_Signal,t);
 fig=figure;
 hold on
 set(fig,'color','white');
 plot(freq,20*log10(abs(Spectrum)),'b','linewidth',2);
 xlabel('Frequency [Hz]');
 ylabel('QPSK Spectrum');
 set(gca,'fontsize',fontsize);
 grid on;
 hold off

%%
%**************** Target description

v=600; %Target speed in meters/second
lambda=c./freq; %Wavelength in meters
dopplershift=(2*v)./lambda; %Convert speed to doppler shift in Hz

doppler_freq=freq+dopplershift; % Create Surveillance signal in Frequency domain
[tempo_surveillance,Surveillance_Signal]=freq2time(Spectrum,doppler_freq);

%**************** Select Surveillance Signals Samples Time Domain 

% Frequency
plot(abs(doppler_freq),abs(Spectrum));
hold on;
plot(abs(freq),abs(Spectrum));
legend('Surveillance Signal','Reference Signal');
xlabel('freq (Hz)');
ylabel('Spectrum');
title('Frequency domain');


%Time
plot(tempo_surveillance,abs(Surveillance_Signal));
hold on;
plot(t,abs(Reference_Signal));
legend('Surveillance Signal','Reference Signal');
xlabel('time (second)');
ylabel('amplitude');
title('Time domain');


k=find(doppler_freq>2.5e7 & doppler_freq<3e7);
Surveillance_Signalcut=Surveillance_Signal(21990:23979);

k2=find(freq>2.5e7 & freq<3e7);
Reference_Signalcut=Spectrum(21990:23979);
%**************** Add White noise to Surveillance_Signal

SNR=20;
Surveillance_Signal=awgn(Surveillance_Signalcut,SNR,'measured');

%**************** Calculate ambiguity and cross-ambiguity functions 

%Reference_Signal ambiguity function
[afmag,doppler] = ambgfun(Reference_Signalcut,fs,1e6,'Cut','Delay');
afmag = afmag*1; % Select plot gain *1
afmag(afmag>1 )= 1;

 %Surveillance_Signal ambiguity function
[afmag2,doppler2] = ambgfun(Surveillance_Signalcut,fs,1e6,'Cut','Delay');
afmag2 = afmag2*1; % Select plot gain *1
afmag2(afmag2>1 )= 1;

%Cross-ambiguity
[afmag3,doppler3] = ambgfun(Reference_Signalcut,Surveillance_Signalcut,fs,[1e6 1e6],'Cut','Delay');
afmag3 = afmag3*1; % Select plot gain *1
afmag3(afmag3>1 )= 1;










