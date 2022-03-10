
% Author of the time2freq function:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the QPSK modulator:
%     Joao Felicio

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022


% The current program is a passive radar simulator, using QPSK.
% For a passive radar: Reference signal=TxQPSK_signal
%                      Surveillance signal=TxQPSK_signal delayed in Time 



% QPSK modulator
clear
close all
clc

%***************** INPUTS ******************
fc=2.45e9; %Carrier frequency
Rb=fc/100; %Bitrate
%SNR=10; %SNR [dB]
fontsize=16;
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

Nb=length(sequence); %Number of bits
Ns=Nb/2; %Number of symbols

%Carrier
Ac=1; %Amplitude
wc=2*pi*fc;
Tc=1/fc; %Period
Tb=1/Rb; %Bit period
fs=4*fc; %Sampling frequency
t_bit=linspace(0,Tb,round(Tb*fs)); %Time base of bit

TxQPSK_signal=[];
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
    TxQPSK_signal=[TxQPSK_signal,QPSK_temp];
end

t=linspace(0,Ns*Tb,length(TxQPSK_signal));


%**************** Plot QPSK signal over time
%fig=figure;
%set(fig,'color','white');
%plot(t,TxQPSK_signal,'linewidth',2,'color','b');
%xlabel('Time [s]');
%ylabel('QPSK signal');
%set(gca,'fontsize',fontsize);
%set(gca,'fontsize',fs);
%grid on;

%**************** Plot frequency spectrum of QPSK signal

 [freq,Spectrum]=time2freq(TxQPSK_signal,t);
 %fig=figure;
 %hold on
 %set(fig,'color','white');
 %plot(freq,20*log10(abs(Spectrum)),'b','linewidth',2);
 %xlabel('Frequency [Hz]');
 %ylabel('QPSK Spectrum');
 %set(gca,'fontsize',fontsize);
 %grid on;
 %hold off



