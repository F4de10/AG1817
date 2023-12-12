% AG1817 HT-22 Rut Hagwall
% ÖVNING 11

clear all
clc
format long g

coordinates_import = importdata("XYZ_20points.txt");
coordinates = coordinates_import.data;

S1x = coordinates(1:20,2);
S1y = coordinates(1:20,3);
S1z = coordinates(1:20,4);

S2x = coordinates(21:40,2);
S2y = coordinates(21:40,3);
S2z = coordinates(21:40,4);

A = [];
L = [];
for i = 1:length(S1x)
    Lx(i) = [S2x(i) - S1x(i)];
    Ly(i) = [S2y(i) - S1y(i)];
    Lz(i) = [S2z(i) - S1z(i)];
    L_temp = [Lx(i); Ly(i); Lz(i)];
    L = [L; L_temp];

    A_temp = [1, 0, 0, S1x(i), 0, -S1z(i), S1y(i);
            0, 1, 0, S1y(i), S1z(i), 0, -S1x(i);
            0, 0, 1, S1z(i), -S1y(i), S1x(i), 0];
    A = [A; A_temp];
    
end

C = eye(60,60);
dX = inv(A' * inv(C) * A) * A' * inv(C) * L;

ds = dX(4) * 10^6;
a1 = dX(5) * (180/pi) * 3600;
a2 = dX(6) * (180/pi) * 3600;
a3 = dX(7) * (180/pi) * 3600;

epsilon = L - A*dX;

n = 20;
sigma_square = (epsilon' * inv(C) * epsilon)/(3*n - 7) % redovisa

Cxx = sigma_square * inv(A' * inv(C) * A);

% redovisa sigma1-7
sigma1 = sqrt(Cxx(1,1))
sigma2 = sqrt(Cxx(2,2))
sigma3 = sqrt(Cxx(3,3))
sigma4 = sqrt(Cxx(4,4)) * 10^6
sigma5 = sqrt(Cxx(5,5)) * (180/pi) * 3600
sigma6 = sqrt(Cxx(6,6)) * (180/pi) * 3600
sigma7 = sqrt(Cxx(7,7)) * (180/pi) * 3600

parametrar = [dX(1); dX(2); dX(3); ds; a1; a2; a3] %redovisa

