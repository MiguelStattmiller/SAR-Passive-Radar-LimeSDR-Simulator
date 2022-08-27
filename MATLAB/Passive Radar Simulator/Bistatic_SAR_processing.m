
% Author of the current program:
% Miguel Albuquerque, Escola Naval, 2022

%The current program includes the steps for bistatic SAR processing, in
%particular:
%. Range compression
%. Range Cell Migration Correction
%. Azimuth Compression

% Inputs: surv_matrix(surveillance signal in frequency domain)
 %        ref_matrix(reference signal in frequency domain)
          



%***************** Range Compression ******************

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


%***************** Curvatures Along Range ******************

[max_correlation, row_idx] = max(range_compressed_matrix, [], 1);
y_max = time_compression_cut(row_idx);
y_max = y_max*c;


%***************** Range Cell Migration Correction ******************


% Calculate necessary time delay between columns, related to the column that has the minimum
% distance


[rows,columns]=size(y_max);
distance_correction=zeros(rows,columns);
for i = 1:size(y_max,1)
thisrow=y_max(i,:);
[R0,ixR0]=min(thisrow(thisrow ~= 0));
first_idx = find(thisrow, 1, 'first');
last_idx  = find(thisrow, 1, 'last');
[rowOfMin, colOfMin] = find(thisrow == R0); 
for col = 1 : columns
        thisColumn = y_max(i, col);
        if thisColumn == 0
            idx2=0;
            continue
        end
        idx2=col-colOfMin;
        if idx2 == 0;
            dist=0;
            continue
        end
        dist=(thisColumn-R0);
        distance_correction(i,col)=dist./1;

end
end

% Distance correction for each column

RMC=y_max-distance_correction;
for i = 1:size(RMC,1)
        thisrow=RMC(i,:);
        isR=find(thisrow);
        RMC2(i,:)=RMC(i,isR);
end


% Range Cell Migration Correction application for each colummn of
% range compressed data.

distance_correction=round(distance_correction);
[rows2,columns2]=size(range_compressed_matrix);
code=0;
    for col2=1:columns2
        thiscol=range_compressed_matrix(:,col2);
        thiscol2=distance_correction(:,col2);
        n=thiscol2;
        migration=circshift(thiscol,n);
        if n>0
            migration(1:n) = 0;
           else
            migration(end+n+1:end) = 0;
        end
   migration_compressed_matrix(:,col2)=migration;
   
    end


%***************** Azimuth compression ******************

% FFT of each Row of migration corrected data
[rows,columns]=size(migration_compressed_matrix);
for row = 1 :  rows
    thisrow = migration_compressed_matrix(row, :);
    [freq_azimuth,azimuth]=time2freq(thisrow,waypoints);
    azimuth_compressed_matrix(row,:)=azimuth;
end

% Select only positive frequencies
k=find(freq_azimuth>0);
freq_azimuth2=freq_azimuth(k);
azimuth_compressed_matrix2=azimuth_compressed_matrix(:,k);

% Interpolation
[rows,columns]=size(azimuth_compressed_matrix);
for row=1:rows
    azimuth2(row,:)=interp1(freq_azimuth2,abs(azimuth_compressed_matrix2(row,:)),linspace(min(freq_azimuth2), ...
    max(freq_azimuth2),length(waypoints)));
end

% Shift of sinc function
[rows,columns]=size(azimuth2);
[maxValue, linearIndexesOfMaxes] = max(range_compressed_matrix(:));
[rowsOfMaxes colsOfMaxes] = find(range_compressed_matrix == maxValue);
[maxValue2, linearIndexesOfMaxes2] = max(azimuth2(:));
[rowsOfMaxes2 colsOfMaxes2] = find(azimuth2 == maxValue2);
shift=abs(colsOfMaxes-colsOfMaxes2);
for row=1:rows
    thisrow=azimuth2(row,:);
 azimuth3(row,:)=circshift(thisrow,shift);
end


% Expected Azimuth Compression

% FFT of each Row of RCMC data
[rows,columns]=size(RMC);
for row = 1 :  rows
    thisrow = RMC(row, :);
    [freq_azimuth,azimuth]=time2freq(thisrow,waypoints);
    azimuth_compressed_matrix(row,:)=azimuth.';
end

%Interpolation
azimuth2=interp1(freq_azimuth,abs(azimuth),linspace(min(freq_azimuth), ...
max(freq_azimuth),length(waypoints)));


%Shift of sinc Function
[min,posin] = find(RMC, 1, 'first');
[min2,posout]  = find(RMC, 1, 'last');
maxvalue3=max(y_max);
[rowsOfMaxes3 colsOfMaxes3] = find(y_max == maxvalue3);
[val,pos] = max(azimuth2);
shift=abs(colsOfMaxes3-pos);
azimuth2=circshift(azimuth2,shift);

[M,I]=max(azimuth2);
y=RMC(I);
x=I;
