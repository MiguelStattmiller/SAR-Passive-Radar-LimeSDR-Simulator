function [P] = LoS(AoI,Lp)
%This function allows to determine if there is line of sight between
%transmitter and a defined pixel
% If there is line of sight, P returns the coordinates of that target
% If there is not line of sight, P=1

% Pixel and Transmitter coordinates without corrleation to Lp

X_transmitter=8; 
Y_transmitter=5;

x_Pixel=1; 
y_Pixel=1;


% Define vector between transmitter and pixel
V=[(X_transmitter-x_Pixel),(Y_transmitter-y_Pixel)];

Npixels=max(abs(V))+1;
coord=0:1/(Npixels-1):1;
reta=zeros(Npixels,2);
reta(:,1)=V(1)*coord+x_Pixel;
reta(:,2)=V(2)*coord+y_Pixel;


 for n=1:Npixels-1
     if AoI(round(reta(n+1,1)),round(reta(n+1,2))) ~= 0
         P=[Lp*round(reta(n+1,1)),Lp*round(reta(n+1,2))]% No line of sight with target coordinates
                break
     else
                 P=1; % Line of sight
                 
            end
            end



 end
     




