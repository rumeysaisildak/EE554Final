function NOMA_THEO(R1_Setup,R2_Setup,pe_Setup,a1_setup,a2_setup)
ps_dB = 0:5:60;
ps = 10.^(0.1*ps_dB);
pr=ps;
%pe_dB=-5;
%pe=10.^(0.1*pe_dB); % Eaves dropper P = Pmax_dB 

ysr=1;
yrd1=1;
yrd2=1;
ysd1=1;
ysd2=1;
yre=1;






for i=1:3
R1=R1_Setup(i);
R2=R2_Setup(i);
%R2=10.^(0.1*R2);
if nargin < 4
a1=0.85;
a2=0.15;

else 
a1= a1_setup(i);
a2= a2_setup(i);    
    
end
pe=10.^(0.1*pe_Setup(i));

%R1=10.^(0.1*R1_db);

Cth1 = 2^(2*R1);
Cth2 = 2^(2*R2);
w = Cth2*pe*a2;
X = Cth2-1;
Elp1=(1-(a2*Cth1))/(a1*a2*Cth1);


A = (1./(a2.*ps.*ysr)) + (1./(a2.*pr.*yrd2));
B = 1./(a2*ps.*ysd2);
gama=((Cth2)).*pe.*yre.*((Cth2))-1 ;%th 1 mi 2 mi belli deðil makalede.



Expression_a1 = (1./((yre.*B.*w)+1)).*exp(-B.*X);
Expression_a2= (1-(exp(-(((B.*w)+(1./(yre)))).*(Elp1./pe))));

Expression_b1 = (1./((yre.*A.*w)+1)).*exp(-A.*X);
Expression_b2=(1-(exp(-((A.*w)+(1/yre)).*(Elp1./pe))));



Expression_c1 = (1./((yre.*w.*A)+(yre.*w.*B)+1)).*exp(-X.*(A+B)); %Son terim belki kuvvet olacak
Expression_c2 = (1-exp(-((w.*A)+(w.*B)+(1./(yre))).*(Elp1./pe)));

Sonuc(i,:)=1-((Expression_a1.*Expression_a2)+(Expression_b1.*Expression_b2)-(Expression_c1.*Expression_c2))


end

for kk=1:3
R1=R1_Setup(kk);
R2=R2_Setup(kk);
R2=10.^(0.1*R2);
pe=10.^(0.1*pe_Setup(kk));
for i=1:length(ps_dB)
Cth1 = 2^(2*R1);
Cth2 = 2^(2*R2);
Elp1=(1-(a2*Cth1))/(a1*a2*Cth1);    
Asymptotic(kk,i) = exp(-((Elp1)/(pe))); 
end
end
% figure
hold on
semilogy(ps_dB, Sonuc(1,:), 'Linewidth', 1,'Color','r','LineStyle','-');
semilogy(ps_dB, Sonuc(2,:), 'Linewidth', 1,'Color','r','LineStyle','-');
semilogy(ps_dB, Sonuc(3,:), 'Linewidth', 1,'Color','r','LineStyle','-');
semilogy(ps_dB, Asymptotic(1,:), 'Linewidth', 1,'Color','r','LineStyle','--');
semilogy(ps_dB, Asymptotic(2,:), 'Linewidth', 1,'Color','r','LineStyle','--');
semilogy(ps_dB, Asymptotic(3,:), 'Linewidth', 1,'Color','r','LineStyle','--');
% grid on
% xlabel('SNR');
% ylabel('Outage Probability');
% axis([0 60  1e-5 1])
end

