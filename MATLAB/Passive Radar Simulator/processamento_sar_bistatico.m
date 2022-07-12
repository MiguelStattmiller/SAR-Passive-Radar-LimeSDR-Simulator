
% Range compression with reference signal

% Normalized values
maxval=max([max(max(abs(surv_matrix))),max(max(abs(ref_matrix)))]);
ref_matrix=ref_matrix/maxval;
surv_matrix=surv_matrix/maxval;
Range_compression=surv_matrix.*conj(ref_matrix); % Range compression FD domain



% Set each column of range compression to Time Domain
[rows,columns]=size(Range_compression);
for col = 1 : columns
    thisColumn = Range_compression(:, col);
    if max(thisColumn)>0
        deebug=1;
    end
    [time_compression,range_compressed]=freq2time(thisColumn, freq_XQPSK_cut);
    idx1=find(time_compression>=0,1);
    idx2=find(time_compression<=7e-6,1,'last');
    time_compression_cut=time_compression(idx1:idx2);
    range_compressed_cut=range_compressed(idx1:idx2);
    range_compressed_matrix(:,col)=range_compressed_cut;
end

range=time_compression_cut*c;

% Range compressed data representation
figure;
maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
contourf(1:400,range,abs(range_compressed_matrix),'edgecolor','none');
%contourf(20*log10(abs(range_compressed_matrix)),'edgecolor','none');
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');


figure;
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
imagesc(1:400,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-20 maxcolorbar]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');


% Identify all the maximum of range compressed
[max_correlation, row_idx] = max(range_compressed_matrix, [], 1);
y_max = time_compression_cut(row_idx);
y_max = y_max*c;
x_max = waypoints;

% Plot the maximum of range compression

Error=(y_max/distance_matrix)*100;
plot(x_max,y_max);
hold on
plot(x_max,distance_matrix);
title('Range over aeroplane movement');
xlabel('Waypoints (m)');
ylabel('Range to target (m)');
legend('Calculated Distance','Estimated Distance');


% azimuth compression 


% Plot the Signals

figure;
imagesc(doppler_freqSurv,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
title('Range Compressed Data');
xlabel('Doppler (Hz)');
ylabel('Range (m)');


n=1;
u=2;
D=sqrt(1-((u.^2*doppler_freqSurv.^2*c.^2)/4*Vr.^2*fc.^2));
temp = ((R1_matrix-Rm_matrix+n*Lp)/c)*fc;
H_BAC=exp(1*j*2*pi*(temp.*D));




