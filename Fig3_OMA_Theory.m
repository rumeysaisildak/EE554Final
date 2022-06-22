function Fig3_OMA_Theory(R1_Setup,R2_Setup,pe_Setup,a1_setup,a2_setup)
ps_dB = 0:5:60;
ps = 10.^(0.1*ps_dB);
pr=ps;

ysr=1;
yrd1=1;
yrd2=1;
ysd1=1;
ysd2=1;
yre=1;



for i=1:3
R1=R1_Setup(i);
R2=R2_Setup(i);
pe=10.^(0.1*pe_Setup(i));
if nargin < 4
a1=0.85;
a2=0.15;

else 
a1= a1_setup(i);
a2= a2_setup(i);    
    
end


Cth1 = 2^(4*R1);
Cth2 = 2^(4*R2);    
    
E1 = (1./((ps).*(ysd1)))+(1./((ps).*(ysr)))+(1./((ps).*(yrd1)));
E2 = (1./((ps).*(ysd2)))+(1./((ps).*(ysr)))+(1./((ps).*(yrd2)));
J1 = Cth1-1;
J2 = Cth2-1;

F1 = E1*Cth1*pe;
F2 = E2*Cth2*pe;

Sonuc_OMA(i,:)= 1-(((1./((yre.*(F1+F2))+1))).*exp((-E1.*J1)-(E2.*J2)))

end

% figure
semilogy(ps_dB, Sonuc_OMA(1,:), 'Linewidth', 1,'Color','b','LineStyle','-');
semilogy(ps_dB, Sonuc_OMA(1,:), 'Marker','o','Linewidth', 1,'Color','#77AC30','LineStyle','none');
% hold on
semilogy(ps_dB, Sonuc_OMA(2,:), 'Linewidth', 1,'Color','b','LineStyle','-');
semilogy(ps_dB, Sonuc_OMA(2,:), 'Marker','hexagram','Linewidth', 1,'Color','#0072BD','LineStyle','none');
semilogy(ps_dB, Sonuc_OMA(3,:), 'Linewidth', 1,'Color','b','LineStyle','-');
semilogy(ps_dB, Sonuc_OMA(3,:), 'Marker','x','Linewidth', 1,'Color','#7E2F8E','LineStyle','none');


% grid on
% xlabel('SNR');
% ylabel('Outage Probability');
% axis([0 60  1e-4 1])
end

