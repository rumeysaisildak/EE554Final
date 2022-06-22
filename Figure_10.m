
clear all;clc ;close all;
N = 1e6; %iterasyon

ps_dB = -10:5:40;
ps=db2pow(ps_dB);

%pe_dB=0;
%pe=db2pow(pe_dB);

pr=ps;

var_hsr=1;
var_hrd1=1;var_hrd2=1;
var_hsd1=1;var_hsd2=1;
var_hre=1;


hsr  = sqrt(var_hsr/2) *(randn(1,N) + 1i*randn(1,N)); % Channel Gain of hsr
hrd1 = sqrt(var_hrd1/2)*(randn(1,N) + 1i*randn(1,N)); % Channel Gain of hrd1
hrd2 = sqrt(var_hrd2/2)*(randn(1,N) + 1i*randn(1,N)); % Channel Gain of hrd2
hsd1 = sqrt(var_hsd1/2)*(randn(1,N) + 1i*randn(1,N));% Channel Gain of hsd1
hsd2 = sqrt(var_hsd2/2)*(randn(1,N) + 1i*randn(1,N));% Channel Gain of hsd2
hre  = sqrt(var_hre/2) *(randn(1,N) + 1i*randn(1,N));% Channel Gain of hre

hsr_square = abs(hsr).^2; % Channel Coefficient hsr
hrd1_square = abs(hrd1).^2; % Channel Coefficient hrd1
hrd2_square = abs(hrd2).^2;% Channel Coefficient hrd2
hsd1_square = abs(hsd1).^2;% Channel Coefficient hsd1
hsd2_square = abs(hsd2).^2;% Channel Coefficient hsd2
hre_square = abs(hre).^2;% Channel Coefficient hre




a1= 0.85;
a2=1-0.85;



pe_Setup = [0 5 10];
for kk = 1:max(size(pe_Setup))

pe=10.^(0.1*pe_Setup(kk));

for i = 1:length(ps_dB)

YSR1_NOMA = (a1*ps(i)*hsr_square)/((a2*ps(i)*hsr_square)+1);
YSR2_NOMA = (a2*ps(i)*hsr_square);

YSD1_NOMA = (a1*ps(i)*hsd1_square)/((a2*ps(i)*hsd1_square)+1);
YSD2_NOMA = (a2*ps(i)*hsd2_square);

YRD1_NOMA = (a1*ps(i)*hrd1_square)/((a2*ps(i)*hrd1_square)+1);
YRD2_NOMA = (a2*ps(i)*hrd2_square);

Tot_Y1 = max(YSD1_NOMA,min(YSR1_NOMA,YRD1_NOMA));
Tot_Y2 = max(YSD2_NOMA,min(YSR2_NOMA,YRD2_NOMA));

YRE1 = (a1*pe*hre_square);
YRE2 = (a2*pe*hre_square);
    
Secrecy_1 = (1+Tot_Y1)./(1+YRE1);
Secrecy_2 = (1+Tot_Y2)./(1+YRE2);

SOP_1 = find(Secrecy_1>1);
SOP_2 = find(Secrecy_2>1);
     
PSOP_1(kk,i)=length(SOP_1)/N;  
PSOP_2(kk,i)=length(SOP_2)/N;  
     
Tot_SOP1(kk,i) =  (PSOP_1(kk,i))+ (PSOP_2(kk,i))-((PSOP_2(kk,i))*(PSOP_1(kk,i)));
             
end
end
figure
semilogy(ps_dB,Tot_SOP1(1,:),'Marker','o','Linewidth',2,'MarkerSize',8,'Color','#77AC30','LineStyle','none'),hold on
semilogy(ps_dB,Tot_SOP1(2,:),'Marker','hexagram','Linewidth',2,'MarkerSize',8,'Color','#0072BD','LineStyle', ' none '),
semilogy(ps_dB,Tot_SOP1(3,:),'Marker','x','Linewidth',2,'MarkerSize',8,'Color','#7E2F8E','LineStyle',' none '),
%SPSC_1=[0.10627, 0.281422, 0.583333, 0.859026, 0.972452, 0.995345, 0.998358,0.998689, 0.998724, 0.998727, 0.998727];
%SPSC_2=[0.0357875, 0.106269, 0.281384, 0.575412, 0.800598, 0.867006,0.877219, 0.878402, 0.878526, 0.878539, 0.87854];
%SPSC_3=[0.0115539, 0.0357877, 0.106211, 0.265494, 0.423794, 0.4769, 0.485458,0.486465, 0.486571, 0.486582, 0.486583];

SPSC_1=[0.190808447554000,	0.475798693921000,	0.841383616660000,	0.989238474406000,	0.999759953530000,	0.999990850325000,	0.999999357730000,	0.999999938752000,	0.999999994692000,	1.00000000000000,	1]
SPSC_2=[0.0659633061880000,	0.190495520254000,	0.470910465076000,	0.815065274674000,	0.968264704441000,	0.995987745325000,	0.999547916400000,	0.999951043700000,	0.999994628040000,	0.999999513984000,	1]
SPSC_3=[0.0216196992640000,	0.0658514650200000,	0.188092238590000,	0.444353841946000,	0.739944678799000,	0.922661991100000,	0.986037970735000,	0.998203226925000,	0.999796731025000,	0.999977390864000,	0.999997945660000]

semilogy(ps_dB, SPSC_1, 'Linewidth', 1,'Color','r','LineStyle','-');
semilogy(ps_dB, SPSC_2, 'Linewidth', 1,'Color','r','LineStyle','-');
semilogy(ps_dB, SPSC_3, 'Linewidth', 1,'Color','r','LineStyle','-');
legend('NOMA λ_E = 0 (dB) sim.','NOMA λ_E = 5 (dB) sim.','NOMA λ_E = 10 (dB) sim.','NOMA ana.', 'FontSize',11)
grid on
xlabel('SNR (dB)','FontSize', 15);
ylabel('SPSC','FontSize', 15);
axis([-10 40  1e-3 1])
