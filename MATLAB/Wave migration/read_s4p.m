function [freq2, S11_fase, S22_fase, S33_fase, S44_fase, S13_fase]=read_s4p(input)

input=fopen('CAL53_0.S4P','r');
cond_neof=1;
n=1;
%skip header
for i=1:1:8
    linha=fgetl(input);
end;

%if format of s4p is RI
if(linha(8)=='R') 

while cond_neof
   linha=fgetl(input);
    
    if linha==-1
        cond_neof=0;
    else
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f%f');
        freq2(n)=clinha(1);       
        S14_fase(n)=clinha(8)+j*clinha(9);   
        S11_fase(n)=clinha(2)+j*clinha(3);  
        S13_fase(n)=clinha(6)+j*clinha(7);
        
        linha=fgetl(input);
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f');
        S24_fase(n)=clinha(7)+j*clinha(8);
        S22_fase(n)=clinha(3)+j*clinha(4);  
        
        linha=fgetl(input);
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f');
        S34_fase(n)=clinha(7)+j*clinha(8);
        S33_fase(n)=clinha(5)+j*clinha(6); 
         
        linha=fgetl(input);
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f');
        S44_fase(n)=clinha(7)+j*clinha(8);
        
        n=n+1;
    end
end
%if format of s4p is dB-angle
else
    while cond_neof
    linha=fgetl(input);
    
    if linha==-1
        cond_neof=0;
    else
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f%f');
        freq2(n)=clinha(1);
        S14_fase(n)=10^(clinha(8)/20)*exp(j*clinha(9)*pi/180);%S23 S24        
        S11_fase(n)=10^(clinha(2)/20)*exp(j*clinha(3)*pi/180);
        S13_fase(n)=10^(clinha(6)/20)*exp(j*clinha(7)*pi/180);
        
        linha=fgetl(input);
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f');
        S24_fase(n)=10^(clinha(7)/20)*exp(j*clinha(8)*pi/180);
        S22_fase(n)=10^(clinha(3)/20)*exp(j*clinha(4)*pi/180); 
        
        linha=fgetl(input);
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f');
        S34_fase(n)=10^(clinha(7)/20)*exp(j*clinha(8)*pi/180);
        S33_fase(n)=10^(clinha(5)/20)*exp(j*clinha(6)*pi/180);
        
        linha=fgetl(input);
        clinha=sscanf(linha,'%f%f%f%f%f%f%f%f');
        S44_fase(n)=10^(clinha(7)/20)*exp(j*clinha(8)*pi/180);
     
        n=n+1;
    end
    end
end
    
fclose(input);


