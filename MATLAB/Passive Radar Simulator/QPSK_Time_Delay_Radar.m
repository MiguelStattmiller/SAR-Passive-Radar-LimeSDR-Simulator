% Author of the time2freq and freq2time functions:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the QPSK modulator:
%     Joao Felicio

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022

% The current program is a passive radar simulator for static targets with zero-doppler values, using QPSK transmitted signal.
% QPSK signal is created based on a message written in .txt file
% For the passive radar: Reference signal= QPSK signal
%                        Surveillance signal=QPSK_signal delayed in Time 



%**************** QPSK Modulator 
clear
close all
clc

%***************** INPUTS QPSK Modulator ******************
fc=30e6; %Carrier frequency
c=3e8; % Speed of light 
lambda=c/fc; %Wavelength
Rb=fc/100; %Bitrate
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

QPSK_Signal=[];
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
    QPSK_Signal=[QPSK_Signal,QPSK_temp];
end

t=linspace(0,Ns*Tb,length(QPSK_Signal));


%**************** Plot QPSK signal over time
fig=figure;
set(fig,'color','white');
plot(t,QPSK_Signal,'linewidth',2,'color','b');
xlabel('Time [s]');
ylabel('QPSK signal [dB]');
set(gca,'fontsize',fontsize);
%set(gca,'fontsize',fs);
grid on;

%**************** Plot frequency spectrum of QPSK signal

 [freq_XQPSK,X_QPSK]=time2freq(QPSK_Signal,t);
 fig=figure;
 hold on
 set(fig,'color','white');
 plot(abs(freq_XQPSK),20*log10(abs(X_QPSK)),'b','linewidth',2);
 xlim ([20e6 40e6]);
 xlabel('Frequency [Hz]');
 ylabel('QPSK Spectrum [dB]');
 set(gca,'fontsize',fontsize);
 grid on;
 hold off

%******************************************** Passive Radar Simulator

%***************** INPUTS  ******************


% Surveillance antenna

alpha_zero_surv=0.68; % Introduce desired receiver antenna angle in radians
D_surv=11; % Antenna directivity in dBi

% Reference antenna

alpha_zeroSVref=pi; % Introduce desired receiver antenna angle in radians
D_ref=11; % Antenna directivity in dBi

% snell function

teta_flutuations=0.5; % Define the flutuation angle value


%**************** Define targets, Surveillance area and radar positions

% White Noise addition 

SNR=-26; %dB


% Define surveillance area and targets
Lp=1; % Pixel length
Nx= zeros(1,200); % dimension in x, horizontal of surveillance area
Ny= zeros(200,1); % dimension in y, vertical  of surveillance area
AoI=Nx.*Ny; % Surveillance area

% Horizontal targets
AoI(1,12) = 1; % define target, set row 4, from column 7-10 to 1, no correlation to Lp
AoI(2,(5:10)) = 1; % define target, set row 4, from column 7-10 to 1, no correlation to Lp
AoI(5,(5:7)) = 1;
AoI(6,(10:12)) = 1;
normal_ntarget1=[1 0]; % Define a normal vector to the target

% Receiver antenna position

% Surveillance antenna
X_receiver=8;
Y_receiver=16;

% Reference antenna

X_receiverref=8;
Y_receiverref=13;

% Transmitter antenna position

X_transmitter=7;


% Aeroplane movement
number_stops=200;
Vr=250; %In meters/second

distance=length(Ny)*Lp; %In meters
total_time=distance/Vr; %In seconds
time_waypoints=total_time/number_stops;
waypoints=[1:Vr*time_waypoints:distance];
waypoints=round(waypoints);


%%
%***************** Processing code  ******************


% Search for targets in surveillance area


for i=1:numel(waypoints)
    Y_transmitter = waypoints(i);

   for xx=1:200
    for yy=1:200
         if AoI(xx,yy) ~= 0 % Target detection
            X_target= xx;
            Y_target= yy;
            VTransmitter_target=[X_target-X_transmitter Y_target-Y_transmitter]; % Vector transmitter-target
            Vectors_product=dot( VTransmitter_target,normal_ntarget1)/norm(VTransmitter_target)*norm(normal_ntarget1);
            angle_transmitter =180-acosd(Vectors_product);
            VTarget_receiver=[X_receiver-X_target Y_receiver-Y_target]; % Vector Target-receiver
            Vectors_product=dot( VTarget_receiver,normal_ntarget1)/norm(VTarget_receiver)*norm(normal_ntarget1);
            angle_receiver =acosd(Vectors_product);
            status=snell_function(angle_transmitter,angle_receiver,teta_flutuations);
            Pr = LoS_receiver(X_receiver,Y_receiver,X_target,Y_target,AoI);
            Pt = LoS_transmitter(X_transmitter,Y_transmitter,X_target,Y_target,AoI);

           
                 if status ==1 & Pt ==1 & Pr ==1

                        Wi_ref = RXGaussian_ref(X_receiverref,Y_receiverref,X_transmitter,Y_transmitter,alpha_zeroSVref,D_ref);
                        Wi_Surv= RXGaussian(X_receiver,Y_receiver,X_target,Y_target,alpha_zero_surv,D_surv);
                        R1=sqrt( (X_transmitter-X_target).^2 + (Y_transmitter-Y_target).^2); % Distance transmitter-target in meters
                        R2=sqrt( (X_receiver-X_target).^2 + (Y_receiver-Y_target).^2); % Distance Receiver-target in meters
                        Rd=sqrt( (X_receiverref-X_transmitter).^2 + (Y_receiverref-Y_transmitter).^2); % Distance Transmitter-Receiver in meters
                        k0=(2*pi*freq_XQPSK)/c; % Wave index variable in radians
                        Vr=Vr*cos(angle_transmitter);
                        lambda2=c./freq_XQPSK; %Wavelength in meters
                        fB=(2*Vr)./lambda2; %Convert speed to doppler shift in Hz
                        Kd=2*pi*fB/c;
                        Surveillance_SignalFD=Wi_Surv*(1/(R1+R2))*X_QPSK.*exp(-1*j*(k0+Kd)*(R1+R2)); % Surveillance Signal frequency domain
                        Surveillance_SignalFD=awgn(Surveillance_SignalFD,SNR,'measured'); % Introduce white gaussian Noise 
                        Reference_SignalFD=Wi_ref*(1/Rd)*X_QPSK.*exp(-1*j*(k0+Kd)*Rd); % Reference Signal frequency domain
                        fprintf('\n Coordenadas do avião(%d,%d)',X_transmitter,Y_transmitter);
                        fprintf('\n Coordenadas do alvo(%d,%d)',X_target,Y_target);
                        fprintf('\n Distância transmissor-alvo R1 %4.2f metros',R1);
                        fprintf('\n Distância alvo-recetor R2 %4.2f metros',R2);
                        fprintf('\n Distância transmissor-recetor Rd %4.2f metros',Rd);
                     continue
                      
                 end   

            end
             
         end
   end
end



%**************** Calculate Reference and Surveillance Signals in Time Domain  

[time_RS,Reference_Signal]=freq2time(Reference_SignalFD,freq_XQPSK);
[time_SS,Surveillance_Signal]=freq2time(Surveillance_SignalFD,freq_XQPSK);


%**************** Select Reference and Surveillance Signals Samples Time Domain 

% Plot signals in Time Domain
plot(time_RS,abs(Reference_Signal));
hold on
plot(time_SS,abs(Surveillance_Signal));
legend('Reference Signal','Surveillance Signal');


% Frequency
plot(abs(freq_XQPSK),abs(Surveillance_SignalFD));
hold on;
plot(abs(freq_XQPSK),abs(Reference_SignalFD));
legend('Surveillance Signal','Reference Signal');
xlabel('freq (Hz)');
ylabel('Spectrum');
title('Frequency domain');


k=find(time_RS>-2e-4 & time_RS<-1.7e-4);
Reference_SignalCut=Reference_Signal(1:179751);

k2=find(time_SS>-2e-4 & time_SS<-1.7e-4);
Surveillance_SignalCut=Surveillance_Signal(1:179751);

%**************** Calculate ambiguity and cross-ambiguity functions 

%Reference_Signal ambiguity function
[afmag,delay] = ambgfun(Reference_SignalCut,fs,1e6,'Cut','Doppler');
afmag = afmag*1; % Select plot gain *1
afmag(afmag>1 )= 1;

 %Surveillance_Signal ambiguity function
[afmag2,delay2] = ambgfun(Surveillance_SignalCut,fs,1e6,'Cut','Doppler');
afmag2 = afmag2*1; % Select plot gain *1
afmag2(afmag2>1 )= 1;

%Cross-ambiguity
[afmag3,delay3] = ambgfun(Reference_SignalCut,Surveillance_SignalCut,fs,[1e6 1e6],'Cut','Doppler');
afmag3 = afmag3*1; % Select plot gain *1
afmag3(afmag3>1 )= 1;


%**************** Plot ambiguity and cross-ambiguity functions


%Plot Ambiguity Function of Sref
[pks1, index] = max(afmag);
xMax = delay(index);
yMax = pks1;
subplot(3,2,1)
plot(delay,afmag,'LineStyle','-.'); 
hold on
%plot3(delay(pks1,index),afmag(pks1,index), pks1, '^r')
textString = sprintf('(%f, %f)', xMax, yMax);
text(xMax, yMax,textString,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');


%Plot Ambiguity Function of Sr
[pks2,index2] = max(afmag2);
xMax2 = delay2(index2);
yMax2 = pks2;
subplot(3,2,2)   
plot(delay2,afmag2,'LineStyle','-'); 
hold on
%plot3(delay2(pks2,index2), afmag2(pks2,index2), pks2, '^r')
textString2 = sprintf('(%f, %f)', xMax2, yMax2);
text(xMax2, yMax2,textString2,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Delay \tau (us)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sr');


% Plot cross-ambiguity function of Sref and Sr

[pks3,index3] = max(afmag3);
xMax3 = delay3(index3);
yMax3 = pks3;
subplot(3,2,3)
plot(delay3,afmag3,'LineStyle','-'); 
hold on
%plot3(delay3(pks3,index3), afmag3(pks3,index3), pks3, '^r')
textString3 = sprintf('(%f, %f)', xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Cross-ambiguity');


% Plot ambiguity function of Sref, Sr and cross-ambiguity

subplot(3,2,4)
plot(delay,afmag,'LineStyle','-.','Color','b'); % Green Sref
hold on
plot(delay2, afmag2,'LineStyle','-','Color','r'); % Red Sr
plot(delay3, afmag3,'LineStyle','--','Color','g'); % blue cross-ambiguity 
hold off
xlim auto;
grid on; 
colorbar;
xlabel('Delay \tau (s)');
ylabel('Ambiguity Function Magnitude');
title('Sref, Sr and cross-ambiguity');
legend('Sref','Sr','cross-ambiguity');

