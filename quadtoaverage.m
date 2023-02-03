% quadtoaverage.m
Xc = load('Xc.txt');
Yc = load('Yc.txt');

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');
Q5 = load('Q5.txt');
Q6 = load('Q6.txt');
Q7 = load('Q7.txt');
Q8 = load('Q8.txt');

Nx = length(Xc);
Ny = length(Yc);

Q1 = reshape(Q1,Ny,Nx)';
Q2 = reshape(Q2,Ny,Nx)';
Q3 = reshape(Q3,Ny,Nx)';
Q4 = reshape(Q4,Ny,Nx)';
Q5 = reshape(Q5,Ny,Nx)';
Q6 = reshape(Q6,Ny,Nx)';
Q7 = reshape(Q7,Ny,Nx)';
Q8 = reshape(Q8,Ny,Nx)';

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
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2 + Qw.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2 + QB3.^2);

lambda(1) = -0.8611363115940525752239465;
lambda(2) = -0.3399810435848562648026658;
lambda(3) = 0.3399810435848562648026658;
lambda(4) = 0.8611363115940525752239465;

weight(1) = 0.3478548451374538573730639;
weight(2) = 0.6521451548625461426269361;
weight(3) = 0.6521451548625461426269361;
weight(4) = 0.3478548451374538573730639;

Nx = length(Xc)/4;
Ny = length(Yc)/4;

Xcc = zeros(1,Nx);
Ycc = zeros(1,Ny);
Q1c = zeros(Nx,Ny);
Q2c = zeros(Nx,Ny);
Q3c = zeros(Nx,Ny);
Q4c = zeros(Nx,Ny);
Q5c = zeros(Nx,Ny);
Q6c = zeros(Nx,Ny);
Q7c = zeros(Nx,Ny);
Q8c = zeros(Nx,Ny);
QPc = zeros(Nx,Ny);
QBPc = zeros(Nx,Ny);
QMachc = zeros(Nx,Ny);

for i = 1:Nx
    Xcc(i) = sum(Xc(4*(i - 1) + 1:4*i))/4;
end

for j = 1:Ny
    Ycc(j) = sum(Yc(4*(j - 1) + 1:4*j))/4;
end

for i = 1:Nx
    for j = 1:Ny
        for i1 = 1:4
            for j1 = 1:4
                Q1c(i,j) = Q1c(i,j) + weight(i1)*weight(j1)*Q1(4*(i - 1) + i1,4*(j - 1) + j1);
                Q2c(i,j) = Q2c(i,j) + weight(i1)*weight(j1)*Q2(4*(i - 1) + i1,4*(j - 1) + j1);
                Q3c(i,j) = Q3c(i,j) + weight(i1)*weight(j1)*Q3(4*(i - 1) + i1,4*(j - 1) + j1);
                Q4c(i,j) = Q4c(i,j) + weight(i1)*weight(j1)*Q4(4*(i - 1) + i1,4*(j - 1) + j1);
                Q5c(i,j) = Q5c(i,j) + weight(i1)*weight(j1)*Q5(4*(i - 1) + i1,4*(j - 1) + j1);
                Q6c(i,j) = Q6c(i,j) + weight(i1)*weight(j1)*Q6(4*(i - 1) + i1,4*(j - 1) + j1);
                Q7c(i,j) = Q7c(i,j) + weight(i1)*weight(j1)*Q7(4*(i - 1) + i1,4*(j - 1) + j1);
                Q8c(i,j) = Q8c(i,j) + weight(i1)*weight(j1)*Q8(4*(i - 1) + i1,4*(j - 1) + j1);
                QPc(i,j) = QPc(i,j) + weight(i1)*weight(j1)*QP(4*(i - 1) + i1,4*(j - 1) + j1);
                QBPc(i,j) = QBPc(i,j) + weight(i1)*weight(j1)*QBP(4*(i - 1) + i1,4*(j - 1) + j1);
                QMachc(i,j) = QMachc(i,j) + weight(i1)*weight(j1)*QMach(4*(i - 1) + i1,4*(j - 1) + j1);
            end
        end
    end
end

xcc = zeros(Nx,Ny);
ycc = zeros(Nx,Ny);

for j = 1:Ny
    xcc(:,j) = Xcc;
end

for i = 1:Nx
    ycc(i,:) = Ycc';
end

Q1 = Q1c;
Q2 = Q2c;
Q3 = Q3c;
Q4 = Q4c;
Q5 = Q5c;
Q6 = Q6c;
Q7 = Q7c;
Q8 = Q8c;
QP = QPc;
QBP = QBPc;
QMach = QMachc;

xc = xcc;
yc = ycc;

drawCloud
                