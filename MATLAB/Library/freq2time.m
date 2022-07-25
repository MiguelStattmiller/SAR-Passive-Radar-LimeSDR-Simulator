function [t,v]=freq2time(V,Fsig)

%
% freq2time - transforms the given frequency charactheristic into a time 
% domain signal
% 
% [u,t,v]=freq2time(V,freq,fc,tau)
%
% Variables input: V - frequency charactheristic given only for positive
%                       frequencies
%                  freq - the positive frequencies for wich the V is given
% Variables output: t - time (both positive and negative part)
%                   v - output pulse
%
% Andela Zaric
% date of creation: 02/09/2012
% date of latest revision: 07/11/2016 (by Joao Felicio)
%

max_f=10*max(Fsig);
f=linspace(-max_f,max_f,2*max_f*(length(Fsig)-1)/(max(Fsig)-min(Fsig))+1);

V_aux=zeros(1,length(f));

V_aux(find(f>=min(Fsig)):find(f>=min(Fsig))+length(Fsig)-1)=V;
V_aux=V_aux+fliplr(conj(V_aux));
f=f(1:end-1);
V_aux=V_aux(1:end-1);

t=linspace(-0.5/(f(2)-f(1)),0.5/(f(2)-f(1)),length(f)+1);
t=t(1:end-1);
v=fftshift(ifft(ifftshift(V_aux)))*sqrt(2*pi)/(t(2)-t(1));

end
