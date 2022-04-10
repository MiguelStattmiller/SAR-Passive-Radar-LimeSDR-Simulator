function [V,x_Pixel,y_Pixel] = raytracing(X_transmitter,Y_transmitter,Lp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

x_Pixel=1;
y_Pixel=5;
X_vetor=x_Pixel-X_transmitter;
Y_vetor=y_Pixel-Y_transmitter;
coordenadas_Vetor= [X_vetor,Y_vetor];
Norma_Vetor= sqrt( (x_Pixel-X_transmitter).^2 + (y_Pixel-Y_transmitter).^2 );
V=[X_vetor/Norma_Vetor,Y_vetor/Norma_Vetor];




end