function [P] = LoS(AoI,X_transmitter,Y_transmitter)
%This function allows to determine if there is line of sight between
%transmitter and a defined pixel
% If there is line of sight, P returns the coordinates of that target
% If there is not line of sight, P=[-1,-1]

% Pixel coordinates
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
                P=[1,1];
                break
     else
                 P=[-1,-1];

            
            end
            end



 end
     




