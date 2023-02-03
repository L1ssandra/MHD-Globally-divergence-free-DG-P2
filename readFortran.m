% drawtest.m

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

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1 = reshape(Q1,Ny,Nx)';
Q2 = reshape(Q2,Ny,Nx)';
Q3 = reshape(Q3,Ny,Nx)';
Q4 = reshape(Q4,Ny,Nx)';
Q5 = reshape(Q5,Ny,Nx)';
Q6 = reshape(Q6,Ny,Nx)';
Q7 = reshape(Q7,Ny,Nx)';
Q8 = reshape(Q8,Ny,Nx)';

%drawRotor
%drawSmoothVortex
%drawOTV
drawall
%drawdiv