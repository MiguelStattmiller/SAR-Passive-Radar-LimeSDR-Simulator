function [status] = snell_function(angle_transmitter,angle_receiver,teta_flutuations)
% This function presents the verification of snell law,for the angle_transmitter (between
% Transmitter and target) and the angle_receiver(between target and
% receiver)
% if verification of this law occurs , status=1 else status=0.

if angle_transmitter <= angle_receiver + teta_flutuations && angle_transmitter >= angle_receiver - teta_flutuations;

    status=1;

else
    status=0;

end