# SAR Passive Radar Simulator| Passive Radar with SDR | MATLAB | GNU RADIO

**OBJECTIVE:**
 The two main goals of this work are:
 
  Develop a passive SAR simulator in MATLAB. In order to study passive SAR evaluation.
  
  Develop a software for real passive radar,based on ISM spectrum (central frequency of 2.45 GHz) on MATLAB using LimeSDR (https://limemicro.com/products/boards/limesdr/), for detection of real targets(aeroplanes in this case).
 
**WORKING PROCESS ABSTRACT:**
 
In order to achieve the main goal, firstly It was used Gnu Radio and Pothos Flow to produce some receivers of ISM/ FM and passive radar examples using LimeSDR (GNU RADIO folder).

Next, It was produced a SAR passive radar simulator (MATLAB folder).

Finally, It was developped the SAR passive radar processing code, for aeroplanes detection(MATLAB folder).

**ORGANIZATION/WORK DESCRIPTION:** 

 **1.LITERATURE-** In this folder there is a document (.PWP) that contains the abstracts with the most important ideas of Scientific Articles/Documents/books/videos used for the investigation of the theme of the Master's Dissertation, as well as the reference to the specific article.

**2.GNU RADIO-** For the first part of the project, it was necessary to install LimeSDR Drivers and libraries , such as Lime Suite, CubicSDR, GNU radio companion, SDR Console,Pothos Flow,SDRangel,...
  Initially it is necessary to connect LimeSDR to your computer,and be able to recognize It as a device. In order to check if this connection was well made, It can be used Lime suite, where you can see LimeSDR status,callibrate Tx and Rx,...

Next, it was used cubicSDR, SDRangel and SDRconsole to evaluate the spectogram of FM/ISM signals received with the receiver antennas and start gaining some knowledge on this theme.

In order to start programming SDR at a high level end, it was used Gnu radio companion and Pothos Flow:

     *Creating a FM receiver in both programs
     
     *Creating a ISM receiver
     
     *Experiment of Passive Radar
     
     *Passive Radar using the LimeSDR

 **3. MATLAB-** In this folder there is the matlab code used.
 
  LIBRARY- All the functions needed for the main program and SAR passive radar simulator.
  
  Main program- Program for the real tests of passive radar using LimeSDR.
  
     *SDR*- Program to receive samples with LimeSDR USB on two different channels.
   
     *SDR_Rx0*- Program to receive samples with LimeSDR USB in channel Rx0.
  
     *SDR_Rx1*- Program to receive samples with LimeSDR USB in channel Rx1.
     
     *SDR_processing*- Program to correlate the acquired samples with LimeSDR USB in both channels. Allowing the detection of the target.
   
  Passive Radar Simulator- Simulator for SAR passive radar using a QPSK modulator.
  
     *QPSK_Time_Delay*- Allows to study the relation between delay and range response in the cross-ambiguity function,with white noise addition for static targets with zero-doppler values.
     
     *QPSK_Doppler*- Allows to study the relation between doppler and range response in the cross-ambiguity function,with white noise addition for non-static targets with doppler values
     
     *Plots*- The code for all the plots and graphics produced.
     
     *Sar_Passive_Radar*- SAR Passive Radar simulator, combinining Delay and doppler responses from targets. Includes Passive SAR processing.

**REFERENCES:**

The code is based on the work from Jockover.

The code is based on the work from RakhDamir.

The code is based on the work from afonsosenica.

The code is based on the work from YAKOTT.
