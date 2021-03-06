clc
clear all
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

x = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0; 0; 0]

q1 = 0.02;
q2 = 0.01;
q3 = 0.01;
r1 = 1000;
r2 = 1000;


Q = [q1 * eye(3) zero zero;
    zero q2 * eye(3) zero;
    zero zero q3 * eye(3)];
R = [r1 * eye(3) zero;
    zero r2 * eye(3)];




b = mf / a ^ 3 * [cos(w0 * t) * sin(im);
    -cos(im);
    2 * sin(w0 * t) * sin(im)];

B = [zero zero
    inv(J), -inv(J) * CrossMatrix(b);
    -inv(Jw) zero];
s = x(1:3);
S = norm(s);

if S>=1
    s = -s/S^2;
    x(1:3) = s;
end
K = lqr(A, B, Q, R);
x1 = [];
t = 100
for j = 1:t+1
    u = - K * x;
    dx = A * x + B * u;
    x = x + dx*.01;
    
    
    s = x(1:3);
    S = norm(s);

    if S>=1
        s = -s/S^2;
        x(1:3) = s;
    end
    x1 = [x1 x(4)];
    
end
y1 = linspace(0,t,t+1);
plot( y1, x1)


                    
                   
                


