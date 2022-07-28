% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022

%This code allows the processing of SDR acquired samples.

clc
clear all
Fs=10e6;

% Import Data

Reference_signal = load('Sinal rececao3_samples.mat');
Reference_signal=cell2mat(struct2cell(Reference_signal));

%Reference_time=load('Canal câmara_timestamp.mat');



Surveillance_signal=load('Aviao5_samples.mat');
Surveillance_signal=cell2mat(struct2cell(Surveillance_signal));

%Surveillance_time=load('Canal 5_timestamp.mat');


% Time base vectors
RS_time=0:1/Fs:numel(Reference_signal)/Fs;
SS_time=0:1/Fs:numel(Surveillance_signal)/Fs;


%**************** Processing time delay code of the received signals

pulse_size = 1000;
nRef = numel(Reference_signal)/pulse_size;
nSurv = numel(Surveillance_signal)/pulse_size;


% Doppler delay processing Sr only
xMax3=[]; yMax3=[];
r=Reference_signal(1:10000000);
counter=1;
surveillance_reshaped = reshape(Surveillance_signal, 1000000, 1, []);
%reference_reshaped = reshape(Reference_signal, 1000000, 1, []);
sz = size(surveillance_reshaped);
for i = 101
    %r = reference_reshaped(:,1,i);
    s = surveillance_reshaped(:,1,i);
    [afmag3,doppler3] = ambgfun(real(r),real(s),Fs,[250000 250000],'Cut','Delay','CutValue',0.4255358*Fs/10000000);
    afmag3(afmag3>1 )= 1;
        if i==41
        dg=1;
    end
    [pks3,index3] = max(afmag3);
    xMax3(:,i)=doppler3(index3);
    yMax3(:,i)=pks3;
    counter=counter+1
    end


    % Doppler delay processing Sr only
xMax3=[]; yMax3=[];
r=Reference_signal(1:10000000);
counter=1;
surveillance_reshaped = reshape(Surveillance_signal, 1000000, 1, []);
%reference_reshaped = reshape(Reference_signal, 1000000, 1, []);
sz = size(surveillance_reshaped);
for i = 101
    %r = reference_reshaped(:,1,i);
    s = surveillance_reshaped(:,1,i);
[afmag,delay] = ambgfun(r,s,Fs,[250000 250000],'cut','Doppler');
afmag(afmag>1 )= 1;
afmag(afmag>1 )= 1;
        if i==101
        dg=1;
    end
    [pks3,index3] = max(afmag);
    xMax3(:,i)=delay(index3);
    yMax3(:,i)=pks3;
    counter=counter+1
    end



% Doppler delay processing
xMax3=[]; yMax3=[];
counter=1;
surveillance_reshaped = reshape(Surveillance_signal, 200000000, 1, []);
reference_reshaped = reshape(Reference_signal, 200000000, 1, []);
sz = size(surveillance_reshaped);
for i = 1:sz(3)
    r = reference_reshaped(:,1,i);
    s = surveillance_reshaped(:,1,i);
    [afmag3,doppler3] = ambgfun(r,s,Fs,[250000 250000],'Cut','Delay');
    afmag3(afmag3>1 )= 1;
    [pks3,index3] = max(afmag3);
    xMax3(:,i)=doppler3(index3);
    yMax3(:,i)=pks3;
    counter=counter+1
    end

figure;
plot(yMax3);
ylabel('Máximo correlação');
xlabel('Número de iterações');
title('Correlação doppler');

figure;
plot(xMax3);
ylabel('Máximo correlação doppler Hz');
xlabel('Número de iterações');
title('Correlação doppler');
ylim([-2e3 2e3]);


for k = nRef
    Reference_signal_samples =Reference_signal((k-1)*pulse_size+(1:pulse_size));

    for kk =nSurv

   
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

r=Surveillance_signal(1:1000);
%reference_reshaped = reshape(Reference_signal, 1000, 4, []);
surveillance_reshaped = reshape(Surveillance_signal, 1000, 4, []);
counter=1;
sz = size(surveillance_reshaped);

sz1=70000;
sz2=75000;
FigH = figure;
AxesH = axes(FigH);
SurfH = surf(AxesH, [], [], [],'LineStyle','-.');  
shading interp;
axis([-13e-6 -11e-6 -2000 2000]); 
grid on; 
view([140,35]); 
colorbar;
ylabel('Doppler (Hz)');
xlabel('Delay(s)');
title('Ambiguity Function Sref');
for i = sz1:sz2
    %r = reference_reshaped(:,1,i); 
    s = surveillance_reshaped(:,1,i);
    [afmag3,delay3,doppler3] = ambgfun(r,s,Fs,[250000 250000]);
    afmag3(afmag3>1 )= 1;
    counter=counter+1
    set(SurfH, 'XData',delay3, 'YData', doppler3, 'ZData', afmag3);
    exportgraphics(AxesH, sprintf('Image%07d.png', counter));

end
