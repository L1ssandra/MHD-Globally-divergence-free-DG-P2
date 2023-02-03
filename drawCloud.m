% drawCloud.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
Qw = Q4./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
QB3 = Q8;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2 + Qw.^2) - 0.5*(QB1.^2 + QB2.^2 + QB3.^2));


M = max(max(Q1));
m = min(min(Q1));
Q1gray = (Q1 - m)/(M - m);
figure(1)
imshow(Q1gray')
title('Density')

MP = max(max(QP));
mP = min(min(QP));
QPgray = (QP - mP)/(MP - mP);
figure(2)
imshow(QPgray')
title('Pressure')

MB1 = max(max(QB1));
mB1 = min(min(QB1));
QB1gray = (QB1 - mB1)/(MB1 - mB1);
figure(3)
imshow(QB1gray')
title('Bx')

MB2 = max(max(QB2));
mB2 = min(min(QB2));
QB2gray = (QB2 - mB2)/(MB2 - mB2);
figure(4)
imshow(QB2gray')
title('By')