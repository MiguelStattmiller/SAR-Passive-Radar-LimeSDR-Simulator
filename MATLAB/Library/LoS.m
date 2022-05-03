function [P,Pr,Pt] = LoS(X_transmitter,Y_transmitter,X_receiver,Y_receiver,X_target,Y_target,AoI)

[P,Pr,Pt]=deal(NaN);
% Define vector between transmitter and pixel
V=[X_transmitter-X_target Y_transmitter-Y_target];
Npixels=max(abs(V))+1;
coord=0:1/(Npixels-1):1;
reta=zeros(Npixels,2);
reta(:,1)=V(1)*coord+X_target;
reta(:,2)=V(2)*coord+Y_target;

% Define vector between receiver and pixel
V2=[X_receiver-X_target Y_receiver-Y_target];
Npixels2=max(abs(V2))+1;
coord2=0:1/(Npixels2-1):1;
reta2=zeros(Npixels2,2);
reta2(:,1)=V2(1)*coord2+X_target;
reta2(:,2)=V2(2)*coord2+Y_target;

% Line of sight transmitter-target

 for n=1:Npixels-1
     if AoI(round(reta(n+1,1)),round(reta(n+1,2))) ~= 0
           Pt=[round(reta(n+1,1)),round(reta(n+1,2))];% No line of sight with target coordinates
           break

     else
           Pt=1; % Line of sight
     end
 end


% Line of sight receiver-target

 for n=1:Npixels-1
     if AoI(round(reta(n+1,1)),round(reta(n+1,2))) ~= 0
           Pr=[round(reta(n+1,1)),round(reta(n+1,2))];% No line of sight with target coordinates
           break

     else
           Pr=1; % Line of sight
     end
 end


if Pt & Pr ==1
    P=1;


end