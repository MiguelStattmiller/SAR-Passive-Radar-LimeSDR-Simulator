
% Author of the time2freq function:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the QPSK modulator:
%     Joao Felicio

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022


% The current program is a passive radar simulator, using QPSK transmitted signal.
% QPSK signal is created based on a message written in .txt file
% For the passive radar: Reference signal= QPSK signal
%                      Surveillance signal=QPSK_signal delayed in Time 



% QPSK modulator
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


%**************** Calculate Surveillance_Signal 

Signal_zeros=padarray(Reference_Signal,[0 1000],0,'both'); %Adds zeros to reference_signal 
Surveillance_Signal=circshift(Signal_zeros,500); %delay Reference signal in time
atraso = finddelay(Reference_Signal,Surveillance_Signal); % Number of samples of delay




%**************** Read data from channels  

%[n,m]=size(Reference_Signal);
%[n2,m2]=size(Surveillance_Signal);


%for column=1:3000:m
     %samples =(Reference_Signal(:,column:column+999));
       % x = transpose(samples);
       % [afmag,delay,doppler] = ambgfun(x,fs,1e6);
       % afmag = afmag*1;
       % afmag(afmag>1 )= 1;
   % end




%**************** Calculate ambiguity and cross-ambiguity functions 

%Reference_Signal ambiguity function
[afmag,delay] = ambgfun(Reference_Signal,fs,1e6,'Cut','Doppler');
afmag = afmag*1;
afmag(afmag>1 )= 1;

 %Surveillance_Signal ambiguity function
[afmag2,delay2] = ambgfun(Surveillance_Signal,fs,1e6,'Cut','Doppler');
afmag2 = afmag2*1;
afmag2(afmag2>1 )= 1;

%Correlation
[afmag3,delay3] = ambgfun(Reference_Signal,Surveillance_Signal,fs,[1e6 1e6],'Cut','Doppler');
afmag3 = afmag3*1;
afmag3(afmag3>1 )= 1;


%**************** Plot ambiguity and cross-ambiguity functions


%Ambiguity Function Sref
[pks1,locs1] = findpeaks(afmag(:));
[r1,c1] = ind2sub(size(afmag), locs1);
[maxValue1] = max(afmag(:));
subplot(3,2,1)
plot(delay,afmag,'LineStyle','-.'); 
hold on
plot3(delay(r1,c1),afmag(r1,c1), pks1, '^r')
hold off
text(-0.5e-9,0.6,maxValue1,['\leftarrow Máximo = ' num2str(maxValue1)],'Color','b','FontSize',10);
shading interp;
xlim ([-8e-5 8e-5]);
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');


%Ambiguity Function Sr
[pks2,locs2] = findpeaks(afmag2(:));
[r2,c2] = ind2sub(size(afmag2), locs2);
[maxValue2] = max(afmag2(:));
subplot(3,2,2)   
plot(delay2,afmag2,'LineStyle','-'); 
hold on
plot3(delay2(r2,c2), afmag2(r2,c2), pks2, '^r')
hold off
text(-0.5e-9,0.6,maxValue2,['\leftarrow Máximo = ' num2str(maxValue2)],'Color','b','FontSize',10);
shading interp;
xlim([-8e-5 8e-5]);
grid on; 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sr');


% Plot the cross-ambiguity function of Sref and Sr

[pks3,locs3] = findpeaks(afmag3(:));
[r3,c3] = ind2sub(size(afmag3), locs3);
[maxValue3] = max(afmag3(:));
subplot(3,2,3)
plot(delay3,afmag3,'LineStyle','-'); 
hold on
plot3(delay3(r3,c3), afmag3(r3,c3), pks3, '^r')
hold off
text(-0.5e-9,0.6,maxValue3,['\leftarrow Máximo = ' num2str(maxValue3)],'Color','b','FontSize',10);
shading interp;
xlim([-8e-5 8e-5]);
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Cross-correlation');


% Plot of Sref, Sr and cross-ambiguity

subplot(3,2,4)
plot(delay,afmag,'LineStyle','-.','Color','g'); % Green Sref
hold on
plot(delay2, afmag2,'LineStyle','-','Color','r'); % Red Sr
plot(delay3, afmag3,'LineStyle','--','Color','b'); % blue corss-ambiguity 
hold off
xlim ([-2e-5 2e-5]);
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Sref, Sr and cross-ambiguity');
legend('Sref','Sr','Cross-ambiguity');
