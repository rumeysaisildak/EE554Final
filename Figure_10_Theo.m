


clc;clear all;close all;
ps_dB = -10:5:40;
ps = 10.^(0.1*ps_dB);
pr=ps;

ysr=1;yrd1=1;
yrd2=1;ysd1=1;
ysd2=1;yre=1;


a1=0.85;
a2=1-a1;
pe_Setup = [0 5 10];

for i=1:3

pe=10.^(0.1*pe_Setup(i)); 
    
A = (1./(a2.*ps.*ysr)) + (1./(a2.*pr.*yrd2));
B = 1./(a2*ps.*ysd2);
qsi = pe*a2;


SPSC_a1=(1./(yre.*B.*qsi + 1));
SPSC_a2=(1 - exp (-((B.*qsi) + (1./yre))*(1./qsi)));


SPSC_b1=(1./(yre.*A.*qsi + 1));
SPSC_b2=(1 - exp (-((A.*qsi) + (1./yre)).*(1./qsi)));


SPSC_c1=(1./((yre.*(A + B).*qsi) + (1)));
SPSC_c2=(1 - exp (-(((A + B).*qsi) + (1./yre)).*(1./qsi)));


Sonuc_SPSC(i,:) = (SPSC_a1.*SPSC_a2)+(SPSC_b1.*SPSC_b2)-(SPSC_c1.*SPSC_c2);
end


figure
semilogy(ps_dB, Sonuc_SPSC(1,:), 'Linewidth', 1,'Color',[0.0,0.0,1.0],'LineStyle','-');
hold on
semilogy(ps_dB, Sonuc_SPSC(2,:), 'Linewidth', 1,'Color',[0.0,0.0,1.0],'LineStyle','-');
semilogy(ps_dB, Sonuc_SPSC(3,:), 'Linewidth', 1,'Color',[0.0,0.0,1.0],'LineStyle','-');

SPSC_1=[0.10627, 0.281422, 0.583333, 0.859026, 0.972452, 0.995345, 0.998358,0.998689, 0.998724, 0.998727, 0.998727];
SPSC_2=[0.0357875, 0.106269, 0.281384, 0.575412, 0.800598, 0.867006,0.877219, 0.878402, 0.878526, 0.878539, 0.87854];
SPSC_3=[0.0115539, 0.0357877, 0.106211, 0.265494, 0.423794, 0.4769, 0.485458,0.486465, 0.486571, 0.486582, 0.486583];
semilogy(ps_dB, SPSC_1, 'Linewidth', 1,'Color',[0.0,0.0,1.0],'LineStyle','-');
semilogy(ps_dB, SPSC_2, 'Linewidth', 1,'Color',[0.0,0.0,1.0],'LineStyle','-');
semilogy(ps_dB, SPSC_3, 'Linewidth', 1,'Color',[0.0,0.0,1.0],'LineStyle','-');

grid on





