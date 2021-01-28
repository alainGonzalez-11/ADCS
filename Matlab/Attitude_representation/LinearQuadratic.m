clc
mu = 3.986005e14;

J1 = 4/120;
J2 = 4/120;
J3 = 8/1200;
Jw1 = 2.65258238E-7;
Jw2 = 2.65258238E-7;
Jw3 = 2.65258238E-7;

al = 657*1000;
im = ((al/1000 + 2.3850e+04) / (57.3 * (900 - 500) / (99 - 97.4)));
a = 6371000+al;
T = 2 * pi * sqrt(a^3/mu)
w0 = 2 * pi / T



J = [J1 0 0; 0 J2 0; 0 0 J3];
Jw = [Jw1 0 0; 0 Jw2 0; 0 0 Jw3];


n = 100;
mf = 7.9e15;

W_RSW = [0, -w0, 0];
zero = [0 0 0; 0 0 0; 0 0 0];
TEMP_1 = CrossMatrix(W_RSW)
TEMP_2 = (4 * CrossMatrix(W_RSW) + 16 * eye(3))
A11 = -0.5*CrossMatrix(W_RSW)  * (4 * CrossMatrix(W_RSW) + 16 * eye(3));

A12 = .25 * eye(3);

A21 = [16 * w0 ^ 2 * (J3 - J2)/J1, 0, 0;
    0, 12 * w0 ^ 2 *(J3 - J1)/J2, 0;
    0, 0, 4 * w0 ^ 2 * (J1 - J2)/J3];

A22 = [0 0 w0 * (J1 - J2 + J3) / J1;
    0 0 0;
    -w0 * (J1 - J2 + J3) / J3 0 0];

A23 = [0 0 w0 * Jw3 /J1;
    0 0 0;
    -w0 * Jw1 / J3 0 0];

A = [A11 A12 zero;
    A21 A22 A23;
    zero zero zero];
t=0
b = (mf / a ^ 3) * [cos(w0 * t) * sin(im);
                        -cos(im);
                        2 * sin(w0 * t) * sin(im)];
B = [zero zero
                        inv(J), -inv(J) * CrossMatrix(b);
                        -inv(Jw) zero];

for j = 0:4
    for j1 = 0:4
        for j2 = 0:4
            for j3 = 0:4
                for j4 = 0:4
                q1 = 0.01* 10^j;
                q2 = 0.01* 10^j1;
                q3 = 0.01* 10^j2;
                r1 = 0.01* 10^j3;
                r2 = 0.01* 10^j4;


                Q = [q1 * eye(3) zero zero;
                    zero q2 * eye(3) zero;
                    zero zero q3 * eye(3)];
                R = [r1 * eye(3) zero;
                    zero r2 * eye(3)];


                Ks = [];


                for i = 1:n

                    t = i * (T/n);
                    b = mf / a ^ 3 * [cos(w0 * t) * sin(im);
                        -cos(im);
                        2 * sin(w0 * t) * sin(im)];

                    B = [zero zero
                        inv(J), -inv(J) * CrossMatrix(b);
                        -inv(Jw) zero];

                    K = lqr(A, B, Q, R);
                    Ks = [Ks; K];
                end
                file = strcat('Ks/', int2str(j), '_', int2str(j1), '_', int2str(j2), '_', int2str(j3), '_', int2str(j4), '.xlsx');
                xlswrite(file, Ks);
                end
            end
        end
    end
end
clear all

