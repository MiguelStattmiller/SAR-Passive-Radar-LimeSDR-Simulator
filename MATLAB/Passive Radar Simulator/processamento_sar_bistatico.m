
% Range compression with reference signal

% Normalized values
maxval=max([max(max(abs(surv_matrix))),max(max(abs(ref_matrix)))]);
ref_matrix=ref_matrix/maxval;
surv_matrix=surv_matrix/maxval;
Range_compression=surv_matrix.*conj(ref_matrix); % Range compression FD domain



% Set each column of range compression to Time Domain
[rows,columns]=size(Range_compression);
maxdistance=600; 
maxtime=maxdistance/c;
range_compressed_matrix=zeros(1600,columns);
for col = 1 : columns
    thisColumn = Range_compression(:, col);
    if max(thisColumn)>0
        deebug=1;
    end
    [time_compression,range_compressed]=freq2time(thisColumn, freq_XQPSK_cut);
    idx1=find(time_compression>=0,1);
    idx2=find(time_compression<=maxtime,1,'last');
    time_compression_cut=time_compression(idx1:idx2);
    range_compressed_cut=range_compressed(idx1:idx2);
    range_compressed_matrix(:,col)=range_compressed_cut.';
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

subplot(3,2,1);
maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
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
maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
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


%RCMC


% Calculate necessary time delay between columns related to the column that has the minimum
% distance

[rows,columns]=size(distance_matrix);
[R0,ixR0]=min(distance_matrix(distance_matrix ~= 0));
first_idx = find(distance_matrix, 1, 'first');
last_idx  = find(distance_matrix, 1, 'last');
[rowOfMin, colOfMin] = find(distance_matrix == R0);  % Find row and col of min.
tshift=zeros(rows,columns);
for col = 1 : columns
        thisColumn = distance_matrix(:, col);
        if thisColumn == 0
            idx2=0;
            continue
        end
        idx2=col-colOfMin;
        if idx2 == 0
            dist=0;
            continue
        end
        dist=thisColumn-R0;
        tshift(:,col)=dist./c;

end

% Distance correction

distance_correction = tshift*c;
RMC=distance_matrix-distance_correction;
expectcol=mean(first_idx:last_idx);
expectvalue=R0;


plot(x_max,RMC);
hold on
plot(x_max,y_max);
%plot(expectcol, expectvalue, '^r');
hold off
title('Range over aeroplane movement');
xlabel('Waypoints (m)');
ylabel('Range to target (m)');
legend('Range Cell Migration correction','Calculated Distance in Range Compression','Expected Distance to Target');

isR=find(y_max);
y_cutted=y_max(isR);
[a,b]=intersect(range,y_cutted);
b=b.';
RMC2=RMC(:,184:255);
range(:,(254:325))=RMC2;


figure;
%subplot(3,2,1);
%maxcolorbar=max(max(20*log10(abs(range_compressed_matrix))));
imagesc(1:400,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-20 maxcolorbar]);
title('Range Migrated Data');
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
%caxis([maxcolorbar-30 maxcolorbar-10]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
ylim([95 125]);

imagesc(azimuth,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar-10]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');


% azimuth compression 


[rows,columns]=size(range_compressed_matrix);
for row = 1 : rows
    thisrow = range_compressed_matrix(row, :);
    [freq_azimuth,azimuth]=time2freq(RMC,waypoints);
end


plot(waypoints,RMC);
hold on
plot(waypoints,azimuth);
hold off
title('Range over aeroplane movement');
xlabel('Waypoints (m)');
ylabel('Range to target (m)');
legend('Range Cell Migration correction','Azimuth compression');
ylim([0 98]);

imagesc(waypoints,range,abs(range_compressed_matrix));
colorbar;
colormap(jet);
%caxis([maxcolorbar-30 maxcolorbar-10]);
title('Range Compressed Data');
xlabel('Waypoints (m)');
ylabel('Range (m)');
