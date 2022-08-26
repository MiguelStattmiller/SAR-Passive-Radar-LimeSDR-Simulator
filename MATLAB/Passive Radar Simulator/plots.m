% The current program represents the code for all the plots needed to comprehend the work developed.


%***************** Spectrum formation of both signals ******************
figure;
ic=find(~all(surv_matrix==0));
hL=plot(abs(doppler_freqSurv),20*log10(abs(surv_matrix(:,ic))));
hLg=legend('1 coluna do avião','2 coluna do avião','3 coluna do avião','4 coluna do avião','5 coluna do avião','6 coluna do avião','7 coluna do avião');
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)');
title('Surveillance Spectrum');


figure;
ic=find(~all(ref_matrix==0));
hL=plot(abs(doppler_freqRef),20*log10(abs(ref_matrix(:,ic))));
hLg=legend('1 coluna do avião','2 coluna do avião','3 coluna do avião','4 coluna do avião','5 coluna do avião','6 coluna do avião','7 coluna do avião');
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)');
title('Reference Spectrum');
%***************** Data Visualization of both signals  ******************

figure;
contour(waypoints,freq_XQPSK_cut,20*log10(abs(surv_matrix)));
colorbar;
xlabel('Waypoints (m)');
ylabel('Frequency (Hz)');
ylim([2.7e7 3.3e7]);
title('Surveillance matrix contour');

figure;
contour(1:400,freq_XQPSK_cut,20*log10(abs(ref_matrix)));
colorbar;
xlabel('Waypoints (m)');
ylabel('frequency (Hz)');
ylim([2.7e7 3.3e7]);
title('Reference matrix contour');

figure;
imagesc(1:400,freq_XQPSK_cut,20*log10(abs(surv_matrix)));
title('Surveillance Signal Raw Data');
colorbar;
xlabel('Waypoints (m)');
ylabel('Frequency Signal (Hz)');
ylim([2.7e7 3.3e7]);
xlim([0 400]);

figure;
imagesc(1:400,freq_XQPSK_cut,20*log10(abs(ref_matrix)));
title('Reference Signal Raw Data');
colorbar;
xlabel('Waypoints (m)');
ylabel('Frequency Signal (Hz)');
ylim([2.7e7 3.3e7]);
xlim([0 400]);

%**************** Plot Reference and Surveillance Signals Time/ Frequency Domain 

% Plot signals in Time Domain
plot(time_RS,abs(Reference_Signal));
hold on
plot(time_SS,abs(Surveillance_Signal));
legend('Reference Signal','Surveillance Signal');
xlabel('Time (s)');
ylabel('Spectrum');
title('Time domain');
xlim([-2e-4 -1.5e-4]);


% Frequency

plot(abs(doppler_freqSurv),abs(Surveillance_SignalFD));
hold on;
plot(abs(doppler_freqRef),abs(Reference_SignalFD));
plot(abs(freq_XQPSK_cut),abs(X_QPSK_cut));
legend('Surveillance Signal','Reference Signal','QPSK');
xlabel('freq (Hz)');
ylabel('Spectrum');
title('Frequency domain');
ylim([0 1.1]);

%***************** Plot the interpolation data  ******************


% Time Domain
plot(RSInterp_time,abs(Reference_interp));
hold on
plot(SSInterp_time,abs(Surveillance_interp));
legend('Reference Signal','Surveillance Signal');
xlabel('Time (s)');
ylabel('Spectrum');
title('Interp Time domain');

%**************** Plot ambiguity and cross-ambiguity functions time delay


%Plot Ambiguity Function of Sref
[pks1, index] = max(afmag);
xMax = delay(index);
yMax = pks1;
subplot(3,2,1)
plot(delay,afmag,'LineStyle','-.'); 
hold on
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
xMax3 =delay3(index3);
yMax3 = pks3;
subplot(3,2,3)
plot(delay3,afmag3,'LineStyle','-'); 
hold on
textString3 = sprintf('(%.2e,%f)',xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
xlim auto;
grid on; 
colorbar;
xlabel('delay (s)');
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
textString3 = sprintf('(%.2e, %f)', xMax3, yMax3);
text(xMax3, yMax3,textString3,"Color",'b','FontSize',10);
hold off
shading interp;
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
%axis([-0.5e-5 0.5e-5 -10 10]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('delay (s)');
ylabel('doppler (Hz)');
zlabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sref');


%Plot Ambiguity Function of Sr
subplot(3,2,2)
surf(delay2,doppler2,afmag2,'LineStyle','-.');
shading interp;
%axis([-0.5e-5 0.5e-5 -10000 10000]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('delay (s)');
ylabel('doppler (Hz)');
zlabel('Ambiguity Function Magnitude');
title('Ambiguity Function Sr');


% Plot cross-ambiguity function of Sref and Sr

figure;
surf(delay3,doppler3,afmag3,'LineStyle','-.');
shading interp;
%axis([-0.5e-5 0.5e-5 -4e6 4e6]); 
grid on; 
view([140,35]); 
colorbar;
xlabel('delay (s)');
ylabel('doppler (Hz)');
zlabel('Ambiguity Function Magnitude');
title('Cross-ambiguity Function');


%***************** Range compressed data representation  ******************
 
figure;
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
contourf(1:400,range,abs(range_compressed_matrix),'edgecolor','none');
%contourf(20*log10(abs(range_compressed_matrix)),'edgecolor','none');
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar]);
title('Range Compressed Data','FontSize',14);
xlabel('Azimuth [m]','FontSize',14);
ylabel('Range [m]','FontSize',14);

subplot(3,2,1);
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
contourf(1:400,range,abs(range_compressed_matrix),'edgecolor','none');
%contourf(20*log10(abs(range_compressed_matrix)),'edgecolor','none');
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
ylim([0 125]);

subplot(3,2,2);
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
contourf(1:400,range,abs(range_compressed_matrix),'edgecolor','none');
%contourf(20*log10(abs(range_compressed_matrix)),'edgecolor','none');
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
ylim([95 125]);



figure;
%subplot(3,2,1);
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
imagesc(1:400,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-20 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');


subplot(3,2,1);
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
imagesc(1:400,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-20 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
ylim([0 125]);

subplot(3,2,2);
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
imagesc(1:400,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-20 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
ylim([95 125]);

%***************** Expected curvatures Along Range ******************


plot(waypoints,y_max);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
legend('Calculated Distance');


plot(waypoints,y_max);
hold on
plot(waypoints,distance_matrix);
title('Range Compressed Data','FontSize',14);
xlabel('Azimuth [m]','FontSize',14);
ylabel('Range [m]','FontSize',14);
legend('Calculated Distance','Estimated Distance'),set(legend,'fontsize',14);



%***************** Expected Range Cell Migration Correction******************



plot(waypoints,RMC);
hold on
plot(waypoints,y_max);
hold off
title('Cell Migration Correction Data','FontSize',14);
xlabel('Azimuth [m]','FontSize',14);
ylabel('Range [m]','FontSize',14);
legend('Cell Migration Correction Data','Calculated Distance','fontsize',14);



%***************** Range Cell Migration Correction in Range Compressed Data******************

figure;
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
contourf(1:400,range,abs(migration_compressed_matrix),'edgecolor','none');
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar-10]);
title('Cell Migration Correction Data','FontSize',14);
xlabel('Azimuth [m]','FontSize',14);
ylabel('Range [m]','FontSize',14);

%***************** Azimuth Compression ******************

contourf(waypoints,range,abs(azimuth3),'edgecolor','none');
colorbar;
colormap(jet);
title('Azimuth Compressed Data','FontSize',14);
xlabel('Azimuth [m]','FontSize',14);
ylabel('Range [m]','FontSize',14);

%***************** Expected Azimuth Compression ******************


figure;
plot(waypoints,azimuth2);
hold on
plot(waypoints,RMC);
plot(x,y,'^r');
title('Target detection','FontSize',14);
xlabel('Azimuth [m]','FontSize',14);
ylabel('Range [m]','FontSize',14);
legend('Azimuth Compression','Cell Range Migration Correction',sprintf('(%2.0f,%2.0f) Target',x,y),'fontsize',14);

