% Author of the time2freq and freq2time functions:
%    Andela Zaric  02/09/2012
%    date of latest revision: 07/11/2016 (by Joao Felicio)

% Author of the QPSK modulator:
%     Joao Felicio

% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022

% The current program is a SAR passive radar simulator for static targets with zero-doppler values, using QPSK transmitted signal.
% QPSK signal is created based on a sequence of bits introduced by the user.
% For the SAR passive radar: 
% .Transmitter is on board of a moving platform.
% .There are 2 static antenna receivers.
% .Transmitted signal= Sequence of bits modulated by a QPSK modulator.

% The simulator contains:
% . Doppler information, delay information, white noise addition, antennas radiation patterns influence, snell function between transmitter and receiver,
% . line of sight between transmitter-target and target-receiver and bistatic SAR processing. 



%**************** QPSK Modulator 
clear
close all
clc

%***************** INPUTS QPSK Modulator ******************
fc=30e6; %Carrier frequency, and the intermidiate frequency for calculations
c=3e8; % Speed of light 
lambda=c/fc; %Wavelength
Rb=fc/100; %Bitrate
fontsize=12;

%*******************************************

% Transmitted Signal uses this sequence of bits

sequence=[1,0,0,1,0,0,1,0,0,1,1,1]; % Introduce the sequence of bits to modulate 
Nb=length(sequence); %Number of bits
Ns=Nb/2; %Number of symbols

%Carrier
Ac=1; %Amplitude
wc=2*pi*fc;
Tc=1/fc; %Period
Tb=1/Rb; %Bit period
fs=4*fc; %Sampling frequency
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
grid on;

%**************** Plot frequency spectrum of QPSK signal

 [freq_XQPSK,X_QPSK]=time2freq(QPSK_Signal,t);
 fig=figure;
 hold on
 set(fig,'color','white');
 plot(abs(freq_XQPSK),20*log10(abs(X_QPSK)),'b','linewidth',2);
 %xlim ([20e6 40e6]);
 xlabel('Frequency [Hz]');
 ylabel('QPSK Spectrum [dB]');
 set(gca,'fontsize',fontsize);
 grid on;
 hold off

%********************************************  SAR Passive Radar Simulator

%***************** INPUTS  ******************

% Frequency RF

F_rf=2.45e9; % Desired central frequency for the simulator

% Bandwith

idxx= find(20*log10(abs(X_QPSK))>-140); % Produce a cut in transmitted spectrum of 30dB
freq_XQPSK_cut=freq_XQPSK(idxx);
freq_XQPSK_cut=abs(freq_XQPSK_cut);
X_QPSK_cut=X_QPSK(idxx);
X_QPSK_cut=abs(X_QPSK_cut);

%********** Plot the part of the spectrum of QPSK signal used to transmit

 fig=figure;
 hold on
 set(fig,'color','white');
 plot(freq_XQPSK_cut,20*log10(X_QPSK_cut),'b','linewidth',2);
 xlabel('Frequency [Hz]');
 ylabel('QPSK Spectrum [dB]');
 set(gca,'fontsize',fontsize);
 grid on;
 hold off

 total_BW=max(freq_XQPSK_cut)-min(freq_XQPSK_cut); % Bandwidth of the transmitted signal

% Transmitting antenna

alpha_zero_transmit=0.2; % Introduce desired receiver antenna angle in radians
D_transmit=7; % Desired Antenna directivity in dBi

theta_3dB_transmit=sqrt(4*pi/(10^(D_transmit/10))); % Transmitter antenna HPBW
teta_flutuations_transm=theta_3dB_transmit/2; % Define the flutuation angle value, used for detection between transmitter-target


% Surveillance antenna

alpha_zero_surv=0.2; % Introduce desired receiver antenna angle in radians
D_surv=5; % Desired Antenna directivity in dBi

% Reference antenna

alpha_zeroSVref=0.87; % Introduce desired receiver antenna angle in radians
D_ref=8; % Desired Antenna directivity in dBi

% snell function

theta_3dB_surv=sqrt(4*pi/(10^(D_surv/10))); % Surveillance Antenna HPBW
teta_flutuations_reflec=theta_3dB_surv/2; % Define the flutuation angle value, used for detection between target-surveillance antenna


%**************** Define targets, Surveillance area and radar positions

% White Noise addition 

SNR=30; % Introduce white noise SNR value, applied to surveillance signal


% Define surveillance area and targets
Lp=1; % Pixel length
Nx= zeros(1,400); % dimension in x, horizontal of surveillance area
Ny= zeros(400,1); % dimension in y, vertical  of surveillance area
AoI=Nx.*Ny; % Surveillance area

% Horizontal targets
AoI(20,(250:300)) = 1; % define targets, set row 4, from column 7-10 to 1.
normal_ntarget1=[0 -1]; % Introduce a normal vector to the targets, used for angles calculations

% Receiver antenna position

% Surveillance antenna coordinates
X_receiver=400;
Y_receiver=400;

% Reference antenna coordinates

X_receiverref=400;
Y_receiverref=400;

% Transmitter antenna x coordinate

X_transmitter=400;


% Transmitter antenna y coordinate, by defining platform moovement along y
% axis

number_stops=400; % Desired Number of times for radar transmission
Vr=250; % platform velocity in meters/second


distance=length(Ny)*Lp; % Area dimension in meters
total_time=distance/Vr; % Time spended by the platform to get across the area dimension in seconds
time_waypoints=total_time/number_stops; % Time between transmissions
waypoints=[1:Vr*time_waypoints:distance]; % Y coordinates for radar transmissions
waypoints=round(waypoints);

%%
%***************** Targets detection code  ******************


surv_matrix=zeros(length(X_QPSK_cut),length(waypoints)); % Define surveillance Signal matrix
ref_matrix=zeros(length(X_QPSK_cut),length(waypoints)); % Define Reference Signal matrix


for i=1:numel(waypoints) % For each waypoint
    Y_transmitter = waypoints(i);

   for xx=1:400
    for yy=1:400
         if AoI(xx,yy) ~= 0 % Target detection
            X_target= xx;
            Y_target= yy;
            VTransmitter_target=[X_target-X_transmitter Y_target-Y_transmitter]; % Vector transmitter-target
            Vectors_product=dot( VTransmitter_target,normal_ntarget1)/(norm(VTransmitter_target)*norm(normal_ntarget1));
            angle_transmitter =180-acosd(Vectors_product); % angle between transmitter-target
            VTarget_receiver=[X_target-X_receiver Y_target-Y_receiver]; % Vector Target-receiver
            Vectors_product=dot( VTarget_receiver,normal_ntarget1)/(norm(VTarget_receiver)*norm(normal_ntarget1));
            angle_receiver =acosd(Vectors_product); % angle between target-receiver
            Vectors_product_2= Colin_vectors(VTransmitter_target,VTarget_receiver);
            status=snell_function(angle_transmitter,angle_receiver,teta_flutuations_reflec,teta_flutuations_transm,Vectors_product_2);
            Pr = LoS_receiver(X_receiver,Y_receiver,X_target,Y_target,AoI); % Line of sight between target-receiver
            Pt = LoS_transmitter(X_transmitter,Y_transmitter,X_target,Y_target,AoI); % Line of sight between transmitter-target

           
                 if status ==1 & Pt ==1 & Pr ==1 % SAR Passive detection

                        Wi_ref = RXGaussian_ref(X_receiverref,Y_receiverref,X_transmitter,Y_transmitter,alpha_zeroSVref,D_ref); % Reference antenna radiation pattern
                        Wi_Surv= RXGaussian(X_receiver,Y_receiver,X_target,Y_target,alpha_zero_surv,D_surv); % Surveillance antenna radiation pattern
                        Wi_transmit=RXGaussian_transmit(X_target,Y_target,X_transmitter,Y_transmitter,alpha_zero_transmit,D_transmit); % Transmitter antenna radiation pattern
                        R1=sqrt( (X_transmitter-X_target).^2 + (Y_transmitter-Y_target).^2); % Distance transmitter-target in meters
                        R2=sqrt( (X_receiver-X_target).^2 + (Y_receiver-Y_target).^2); % Distance Receiver-target in meters
                        Rd=sqrt( (X_receiverref-X_transmitter).^2 + (Y_receiverref-Y_transmitter).^2); % Distance Transmitter-Receiver in meters
                        k0=(2*pi*freq_XQPSK_cut)/c; % Wave index variable in radians
                        V=Vr*cosd(90-angle_transmitter); % Platform velocity for doppler information
                        lambda2=c./F_rf; % Wavelength in meters
                        fB=V./lambda2; % Convert speed to doppler shift in Hz
                        Kd=(2*pi*fB)/c; % Surveillance Doppler information
                        fB_ref=Vr./lambda2; % Convert speed to doppler shift in Hz
                        kd_ref=(2*pi*fB_ref)/c; % Reference Doppler information
                        doppler_freqSurv=freq_XQPSK_cut+fB; % Surveillance Doppler frequency
                        doppler_freqRef=freq_XQPSK_cut+fB_ref; % Reference Doppler frequency
                        Surveillance_SignalFD=Wi_Surv.*(1/(R1+R2))*X_QPSK_cut.*exp(-1*j*(k0+Kd)*(R1+R2)); % Surveillance Signal frequency domain
                        Surveillance_SignalFD=awgn(Surveillance_SignalFD,SNR,'measured'); % Introduce white gaussian Noise
                        surv_matrix(:,i)=Surveillance_SignalFD;  % Surveillance Signal of entire detections in frequency domain
                        Reference_SignalFD=Wi_ref.*(1/Rd)*X_QPSK_cut.*exp(-1*j*(k0+kd_ref)*Rd); % Reference Signal frequency domain
                        ref_matrix(:,i)=Reference_SignalFD; % Reference Signal of entire detections in frequency domain
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



%***************** Spectrum formation of both signals ******************
figure;
ic=find(~all(surv_matrix==0));
hL=plot(abs(doppler_freqSurv),abs(surv_matrix(:,ic)));
hLg=legend('1 coluna do avião','2 coluna do avião','3 coluna do avião','4 coluna do avião','5 coluna do avião','6 coluna do avião','7 coluna do avião');
title('Surveillance Spectrum');

figure;
ic=find(~all(ref_matrix==0));
hL=plot(abs(doppler_freqRef),abs(ref_matrix(:,ic)));
hLg=legend('1 coluna do avião','2 coluna do avião','3 coluna do avião','4 coluna do avião','5 coluna do avião','6 coluna do avião','7 coluna do avião');
title('Reference Spectrum');
%***************** Data Visualization of both signals  ******************

figure;
contour(waypoints,doppler_freqSurv,20*log10(abs(surv_matrix)));
colorbar;
xlabel('Waypoints (m)');
ylabel('Frequency (Hz)');
title('Surveillance matrix contour');

figure;
contour(waypoints,doppler_freqRef,20*log10(abs(ref_matrix)));
colorbar;
xlabel('Waypoints (m)');
ylabel('frequency (Hz)');
ylim([2.7e7 3.3e7])
title('Reference matrix contour');

figure;
imagesc(waypoints,doppler_freqSurv,20*log10(abs(surv_matrix)));
title('Surveillance Signal Raw Data');
colorbar;
xlabel('Waypoints (m)');
ylabel('Frequency Signal (Hz)');
xlim([0 400]);

figure;
imagesc(waypoints,doppler_freqRef,20*log10(abs(ref_matrix)));
title('Reference Signal Raw Data');
colorbar;
xlabel('Waypoints (m)');
ylabel('Frequency Signal (Hz)');
xlim([0 400]);

%***************** Signals Time Domain Calculation  ******************

[n, m]   = size(surv_matrix);
[n2, m2] = size(ref_matrix);

 [time_SS,Surveillance_Signal]=freq2time(m,doppler_freqSurv);

 [time_RS,Reference_Signal]=freq2time(m2,doppler_freqRef);




%***************** Signals Correlation Calculation  ******************



afmag3   = cell(m, m2);
delay3   = cell(m, m2);
doppler3 = cell(m, m2);
for i=1:m
    S = surv_matrix(:,i);      % Surveillance of cut signal
    for ii=1:m2
        R = ref_matrix(:,ii);  % Reference of cut signal
        [tmp, delay3{i, ii}, doppler3{i, ii}] = ambgfun(R, S, fs, [1e6 1e6]);
        afmag3{i, ii} = max(1, tmp);
    end
end



%**************** Calculate Reference and Surveillance Signals in Time Domain  

[time_RS,Reference_Signal]=freq2time(Reference_SignalFD,doppler_freqRef);
[time_SS,Surveillance_Signal]=freq2time(Surveillance_SignalFD,doppler_freqSurv);


%**************** Select Reference and Surveillance Signals Samples Time Domain 

% Plot signals in Time Domain
plot(time_RS,abs(Reference_Signal));
hold on
plot(time_SS,abs(Surveillance_Signal));
legend('Reference Signal','Surveillance Signal');


% Frequency
plot(abs(doppler_freqSurv),abs(Surveillance_SignalFD));
hold on;
plot(abs(doppler_freqRef),abs(Reference_SignalFD));
plot(abs(freq_XQPSK_cut),abs(X_QPSK_cut));
legend('Surveillance Signal','Reference Signal','QPSK');
xlabel('freq (Hz)');
ylabel('Spectrum');
title('Frequency domain');

% Time delay cut

k=find(time_RS>-1.95e-4 & time_RS<-1.89e-4);
Reference_SignalCut=Reference_Signal(29869:65838);

k2=find(time_SS>-1.95e-4 & time_SS<-1.89e-4);
Surveillance_SignalCut=Surveillance_Signal(29869:65838);



%************* Calculate ambiguity and cross-ambiguity functions-time delay

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


%*********** Calculate ambiguity and cross-ambiguity functions-doppler delay


%Reference_Signal ambiguity function
[afmag,doppler] = ambgfun(Reference_SignalCut,fs,1e6,'Cut','Delay');
afmag = afmag*1; % Select plot gain *1
afmag(afmag>1 )= 1;

 %Surveillance_Signal ambiguity function
[afmag2,doppler2] = ambgfun(Surveillance_SignalCut,fs,1e6,'Cut','Delay');
afmag2 = afmag2*1; % Select plot gain *1
afmag2(afmag2>1 )= 1;

%Cross-ambiguity
[afmag3,doppler3] = ambgfun(Reference_SignalCut,Surveillance_SignalCut,fs,[1e6 1e6],'Cut','Delay');
afmag3 = afmag3*1; % Select plot gain *1
afmag3(afmag3>1 )= 1;

%************* Calculate ambiguity and cross-ambiguity functions-time and doppler delay

%Reference_Signal ambiguity function
[afmag,delay,doppler] = ambgfun(Reference_SignalCut,fs,1e6);
afmag = afmag*1; % Select plot gain *1
afmag(afmag>1 )= 1;

 %Surveillance_Signal ambiguity function
[afmag2,delay2,doppler2] = ambgfun(Surveillance_SignalCut,fs,1e6);
afmag2 = afmag2*1; % Select plot gain *1
afmag2(afmag2>1 )= 1;

%Cross-ambiguity
[afmag3,delay3,doppler3] = ambgfun(Reference_SignalCut,Surveillance_SignalCut,fs,[1e6 1e6]);
afmag3 = afmag3*1; % Select plot gain *1
afmag3(afmag3>1 )= 1;


%**************** Plot ambiguity and cross-ambiguity functions time delay


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



%**************** Plot ambiguity and cross-ambiguity functions doppler delay

%Plot Ambiguity Function of Sref
[pks1, index] = max(afmag);
xMax = doppler(index);
yMax = pks1;
subplot(3,2,1)
plot(doppler,afmag,'LineStyle','-.'); 
hold on
%plot3(doppler(pks1,index),afmag(pks1,index), pks1, '^r')
textString = sprintf('(%f, %f)', xMax, yMax);
text(xMax, yMax,textString,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');


%Plot Ambiguity Function of Sr
[pks2,index2] = max(afmag2);
xMax2 = doppler2(index2);
yMax2 = pks2;
subplot(3,2,2)   
plot(doppler2,afmag2,'LineStyle','-'); 
hold on
%plot3(doppler2(pks2,index2), afmag2(pks2,index2), pks2, '^r')
textString2 = sprintf('(%f, %f)', xMax2, yMax2);
text(xMax2, yMax2,textString2,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sr');


% Plot cross-ambiguity function of Sref and Sr

[pks3,index3] = max(afmag3);
xMax3 = doppler3(index3);
yMax3 = pks3;
subplot(3,2,3)
plot(doppler3,afmag3,'LineStyle','-'); 
hold on
textString3 = sprintf('(%f, %f)', xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Cross-ambiguity');


% Plot ambiguity function of Sref, Sr and cross-ambiguity

subplot(3,2,4)
plot(doppler,afmag,'LineStyle','-.','Color','b'); % Green Sref
hold on
plot(doppler2, afmag2,'LineStyle','-','Color','r'); % Red Sr
plot(doppler3, afmag3,'LineStyle','--','Color','g'); % blue cross-ambiguity 
hold off
xlim ([0.5e6 1.2e6]);
grid on; 
colorbar;
xlabel('Doppler (Hz)');
ylabel('Ambiguity Function Magnitude');
title('Sref, Sr and cross-ambiguity');
legend('Sref','Sr','cross-ambiguity');

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




