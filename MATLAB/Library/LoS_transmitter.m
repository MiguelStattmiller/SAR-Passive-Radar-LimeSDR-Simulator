function [Pt] = LoS_transmitter(X_transmitter,Y_transmitter,X_target,Y_target,AoI)
%This function allows to determine if there is line of sight between
%transmitter and a target
% If there is not line of sight, P returns the coordinates of that target
% If there is  line of sight, P=1

% Define vector between transmitter and pixel
V=[X_transmitter-X_target Y_transmitter-Y_target];


Npixels=max(abs(V))+1;
coord=0:1/(Npixels-1):1;
reta=zeros(Npixels,2);
reta(:,1)=V(1)*coord+X_target;
reta(:,2)=V(2)*coord+Y_target;



% Line of sight transmitter-target

 for n=1:Npixels-1
     if AoI(round(reta(n+1,1)),round(reta(n+1,2))) ~= 0
           Pt=[round(reta(n+1,1)),round(reta(n+1,2))];% No line of sight with target coordinates
           break

     else
                 Pt=1; % Line of sight
     end
 end

          
    
  


 





