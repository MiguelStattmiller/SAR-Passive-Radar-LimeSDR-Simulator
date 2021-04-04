fc=3.6e9;
c=3e8;
lambda=c/fc;

%Retirar as reflexões internas da antena
%LEITURA PARA (0,0,0)- MEDIDAS 51

Salvo=read_s4p('TESTE53_0.S4P');
Scal=read_s4p('CAL53_0.S4P');
s11=Salvo-Scal;

%Definir a grelha de pesquisa
Nx=41;
Ny=41;
x= -20:20;
y= 0:40;
z= 0:40;
%Posições da antena
ya=0;
za=0;
xa= -50:50;

%Distancia entre antena e pixel
n=3;
xp=[x zeros(size(x,1),60)]; %adiciona zeros para ter a mesma dimensão da matriz x
yp=[y zeros(size(x,1),60)]; %adiciona zeros para ter a mesma dimensão da matriz y
zp=[z zeros(size(x,1),60)]; %adiciona zeros para ter a mesma dimensão da matriz z
dtotal=sqrt( (xa-xp).^2 + (ya-yp).^2 + (za-zp).^2);

%Atraso de fase
phi= (2*pi/lambda)* dtotal;
phip=[phi zeros(size(x,1),1500)]; %adiciona zeros para ter a mesma dimensão da matriz s11
phitotal= s11.*exp(2.* phip .*1i);

for xa= -50:50
    for x= -20:20
        for y= 0:40
            for z= 0:40
            dtotal;
            phitotal;
            image(xa,x,y);
            image(phitotal);
          
        end
    end
  end
end
    



        
