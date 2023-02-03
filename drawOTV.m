% drawOTV.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q5;
QB1 = Q6;
QB2 = Q7;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

figure(1);
contour(xc,yc,Qrho,15);
%mesh(xc,yc,Q1);
%title('rho')
colormap(cool);

figure(2)
contour(xc,yc,QP,15);
%mesh(xc,yc,Q1);
%title('rho')
colormap(cool);
