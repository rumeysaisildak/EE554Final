function NOMA_THEO_Fig7(pe_Setup)
ps_dB = 40;
ps = 10.^(0.1*ps_dB);
pr=ps;

ysr=1;yrd1=1;
yrd2=1;ysd1=1;
ysd2=1;yre=1;

a1=0.15:0.05:0.9;
a2=1-a1;
R1=0.1;
R2=0.5;
%pe_Setup = [0 4 7];

for kk=1:3;
pe=10.^(0.1*pe_Setup(kk));

for i=1:length(a1)
Cth1 = 2^(2*R1);
Cth2 = 2^(2*R2);
w = Cth2*pe*a2(i);
X = Cth2-1;
Elp1=(1-(a2(i)*Cth1))/(a1(i)*a2(i)*Cth1);


A = (1./(a2(i)*ps*ysr)) + (1./(a2(i)*pr*yrd2));
B = 1./(a2(i)*ps.*ysd2);
gama=((Cth2)).*pe.*yre.*((Cth2))-1 ;%th 1 mi 2 mi belli deðil makalede.



Expression_a1 = (1/((yre.*B.*w)+1))*exp(-B.*X);
Expression_a2= (1-(exp(-(((B.*w)+(1/(yre))))*(Elp1/pe))));

Expression_b1 = (1/((yre*A*w)+1))*exp(-A*X);
Expression_b2=(1-(exp(-((A*w)+(1/yre))*(Elp1/pe))));

Expression_c1 = (1./((yre*w*A)+(yre*w.*B)+1))*exp(-X*(A+B)); %Son terim belki kuvvet olacak
Expression_c2 = (1-exp(-((w*A)+(w*B)+(1/(yre)))*(Elp1/pe)));

Sonuc(kk,i)=1-((Expression_a1*Expression_a2)+(Expression_b1*Expression_b2)-(Expression_c1*Expression_c2))


    
end

end
% for i=1:3
% R1=R1_Setup(i);
% R2=R2_Setup(i);
% %R2=10.^(0.1*R2);
% 
% pe=10.^(0.1*pe_Setup(i));
% 
% %R1=10.^(0.1*R1_db);
% 
% Cth1 = 2^(2*R1);
% Cth2 = 2^(2*R2);
% w = Cth2.*pe.*a2;
% X = Cth2-1;
% Elp1=(1-(a2.*Cth1))/(a1.*a2.*Cth1);
% 
% 
% A = (1./(a2.*ps.*ysr)) + (1./(a2.*pr.*yrd2));
% B = 1./(a2.*ps.*ysd2);
% gama=((Cth2)).*pe.*yre.*((Cth2))-1 ;%th 1 mi 2 mi belli deðil makalede.
% 
% 
% 
% Expression_a1 = (1./((yre.*B.*w)+1)).*exp(-B.*X);
% Expression_a2= (1-(exp(-(((B.*w)+(1./(yre)))).*(Elp1./pe))));
% 
% Expression_b1 = (1./((yre.*A.*w)+1)).*exp(-A.*X);
% Expression_b2=(1-(exp(-((A.*w)+(1/yre)).*(Elp1./pe))));
% 
% 
% 
% Expression_c1 = (1./((yre.*w.*A)+(yre.*w.*B)+1)).*exp(-X.*(A+B)); %Son terim belki kuvvet olacak
% Expression_c2 = (1-exp(-((w.*A)+(w.*B)+(1./(yre))).*(Elp1./pe)));
% 
% Sonuc(i,:)=1-((Expression_a1.*Expression_a2)+(Expression_b1.*Expression_b2)-(Expression_c1.*Expression_c2))
% 
% 
% end


hold on
plot(a1, Sonuc(1,:), 'Linewidth', 1,'Color','r','LineStyle','-');
plot(a1, Sonuc(2,:), 'Linewidth', 1,'Color','r','LineStyle','-');
plot(a1, Sonuc(3,:), 'Linewidth', 1,'Color','r','LineStyle','-');
 grid on
% xlabel('SNR');
% ylabel('Outage Probability');
 %axis([0.15 0.9  0 1])

end
