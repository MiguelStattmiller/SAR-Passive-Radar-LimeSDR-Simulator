function [status] = snell_function(angle_transmitter,angle_receiver,teta_flutuations_reflec,teta_flutuations_transm,Vectors_product_2)
% This function presents the verification of snell law,for the angle_transmitter (between
% Transmitter and target) and the angle_receiver(between target and
% receiver)
% if verification of this law occurs , status=1 else status=0.

angle_transmitter=deg2rad(angle_transmitter);
angle_receiver=deg2rad(angle_receiver);

if angle_transmitter - teta_flutuations_transm <= angle_transmitter && angle_transmitter <= angle_transmitter + teta_flutuations_transm &&  Vectors_product_2 ~= 1 && angle_receiver - teta_flutuations_reflec <= angle_receiver  && angle_receiver <= angle_receiver + teta_flutuations_reflec

    status=1;

else 
    status=0;
    
end



