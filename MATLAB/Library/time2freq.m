function [ freq, U ] = time2freq( u_in, t_in )
%
% time2freq - transforms the given time charactheristic into a frequency 
% domain signal
% 
% [freq,U]=time2freq(u, t)
%
% Variables input: u_in - undistorted pulse in TD
%                  t_in - time
% Variables output: U - Gaussian modulated pulse in FD (both positive and negative part)
%                   freq - time (both positive and negative part)
%
fs=1/(t_in(2)-t_in(1));

% t_max=1e-7;
t_max=max(t_in)*10;

t=linspace(-t_max,t_max,2*t_max*fs-1);
u=zeros(1,length(t));

u(find(t>=min(t_in),1):find(t>=max(t_in),1)-1)=interp1(t_in,u_in,t(find(t>=min(t_in),1):find(t>=max(t_in),1)-1));

% freq=linspace(-fs/2,fs/2,length(t));
freq=linspace(-fs/2,fs/2,length(t));

U=zeros(1,length(freq));

U=fftshift(fft(u)*(abs(t(2)-t(1))/sqrt(2*pi)));

end

