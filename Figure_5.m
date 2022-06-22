
clear all;clc ;close all;
N = 1e6; %iterasyon

ps_dB = 0:5:60;
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

%R1,R2,Pe

R1_Setup = [ 0.5 0.5 0.1];
R2_Setup = [  2  2  0.5];
pe_Setup = [  0  -2  -2];
for kk = 1:max(size(R1_Setup))
R1=R1_Setup(kk);
R2=R2_Setup(kk);
%R2=10.^(0.1*R2);
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
    
Capacity_SD1 = 0.5*log2(1+Tot_Y1);
Capacity_SD2 = 0.5*log2(1+Tot_Y2);


Capacity_RE1 = 0.5*log2(1+YRE1);
Capacity_RE2 = 0.5*log2(1+YRE2);

%Secrecy_NOMA1= 0.5*log((1+Tot_Y1)./(1+YRE1));
%Secrecy_NOMA2= 0.5*log((1+Tot_Y2)./(1+YRE2));

Secrecy_NOMA1 = max(0,Capacity_SD1-Capacity_RE1);
Secrecy_NOMA2 =  max(0,Capacity_SD2-Capacity_RE2);


    SOP_1 = find(Secrecy_NOMA1>R1);
    SOP_2 = find(Secrecy_NOMA2>R2);
     
     PSOP_1(kk,i)=length(SOP_1)/N;  
     PSOP_2(kk,i)=length(SOP_2)/N;  
     
     %Tot_SOP(kk,i) = 1-((PSOP_1(kk,i))*PSOP_2(kk,i));
     
     Tot_SOP1(kk,i) =  (1-PSOP_1(kk,i))+ (1-PSOP_2(kk,i))-((1-PSOP_2(kk,i))*(1-PSOP_1(kk,i)));
           
     
end
end
figure
semilogy(ps_dB,Tot_SOP1(1,:),'Marker','o','Linewidth',2,'MarkerSize',8,'Color','#77AC30','LineStyle','none'),hold on
semilogy(ps_dB,Tot_SOP1(2,:),'Marker','hexagram','Linewidth',2,'MarkerSize',8,'Color','#0072BD','LineStyle', ' none '),
semilogy(ps_dB,Tot_SOP1(3,:),'Marker','x','Linewidth',2,'MarkerSize',8,'Color','#7E2F8E','LineStyle',' none '),

NOMA_THEO(R1_Setup,R2_Setup,pe_Setup);
Fig3_OMA_Theory(R1_Setup,R2_Setup,pe_Setup);
legend('R_1=0.5, R_2=2, λ_E=0 (dB) sim.','R_1=0.5, R_2=2, λ_E=-2 (dB) sim.','R_1=0.1, R_2=0.5, λ_E=-2 (dB) sim.', 'NOMA ana.', '', '', 'NOMA asymptotic', '', '', 'OMA ana.', '', '', 'FontSize',11)

grid on
xlabel('SNR (dB)', 'FontSize', 15);
ylabel('Secrecy Outage Probability', 'FontSize', 15);
axis([0 60  1e-4 1])
