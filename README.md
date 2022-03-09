# Passive Radar| SDR | MATLAB | GNU RADIO

OBJECTIVE:

 The main purpose of this work is to develop a passive radar with MATLAB using LimeSDR ( based on ISM spectrum, central frequency of 2.45 GHz).
 SDR used: https://limemicro.com/products/boards/limesdr/
 
 
 WORKING PROCESS:
 
In order to achieve that goal, firstly It was used Gnu Radio and Pothos Flow to produce some receivers of ISM/ FM and passive radar examples using LimeSDR (GNU RADIO folder).
Next, It was developped the main processing code in MATLAB (MATLAB folder).
Finally, It was produced a passive radar simulator (MATLAB folder).


ORGANIZATION AND SEQUENCE: 

 1. LITERATURE- In this folder there is a document (.PWP) that contains the abstracts of Scientific Articles/Documents used for the investigation of the theme of the Master's
Dissertation, as well as the reference to the article used.

 2. GNU RADIO- For the first part of the project, it was necessary to install LimeSDR Drivers and libraries , such as Lime Suite, CubicSDR, GNU radio companion, SDR Console,Pothos Flow and SDRangel.
Initially it is necessary to use Lime suite in order to connect LimeSDR and callibrate Tx and Rx.
Next, it was used cubicSDR, SDRangel, SDRconsole to evaluate the spectogram of FM and ISM signals and start gaining some knowledge on this theme.
In order to start programming SDR at a high level end, it was used Gnu radio companion and Pothos Flow:
     |Creating a FM receiver in both programs
     |Creating a ISM receiver
     |Experiment of Passive Radar
     |Passive Radar using the LimeSDR

 3. MATLAB- In this folder there is the matlab code used.
  LIBRARY- All the functions needed for the main program (All the functions are based on the work from Jockover and RakhDamir).
  PASSIVE RADAR- Main program of passive radar using LimeSDR.


REFERENCES:

The code is based on the work from Jockover.
The code is based on the work from RakhDamir.
The code is based on the work from afonsosenica.
The code is based on the work from YAKOTT.
