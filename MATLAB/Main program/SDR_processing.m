clc
clear all

% Inputs

Fs= 10e6; % Frequency of sampling frequency, Hz   

% Import Data

Reference_signal = load('Canal zero_samples.mat');
Reference_signal=cell2mat(struct2cell(Reference_signal));

Reference_time=load('Canal zero_timestamp.mat');



Surveillance_signal=load('Canal um_samples.mat');
Surveillance_signal=cell2mat(struct2cell(Surveillance_signal));

Surveillance_time=load('Canal um_timestamp.mat');

% Reference Signal Time Base

val=['19-05-21_475'
     '19-05-33_280'];
val = string(val);
val = replace(val, '-', ':');
val = replace(val, '_', '.');
d = duration(val);
d.Format = 'hh:mm:ss.SSS';
RS_time=seconds(diff(d)); % Difference in seconds


% Surveillance Signal Time Base

val2=['17-10-58_086'
    '17-11-09_923'];
val2 = string(val2);
val2 = replace(val2, '-', ':');
val2 = replace(val2, '_', '.');
d2 = duration(val2);
d2.Format = 'hh:mm:ss.SSS';
SS_time=seconds(diff(d2)); % Difference in seconds


% Interpolation

RSInterp_time=[RS_time(1):16/Fs:RS_time(end)];
Reference_interp=interp1(RS_time,Reference_signal,RSInterp_time);
% RSInterp_time time reference for sinal reference_interp 

SSInterp_time=[time_SS(1):16/fs:time_SS(end)];
Surveillance_interp=interp1(time_SS,Surveillance_Signal,SSInterp_time);
% SSInterp_time time reference for sinal Surveillance_interp 

fs=1/(RSInterp_time(2)-RSInterp_time(1));
