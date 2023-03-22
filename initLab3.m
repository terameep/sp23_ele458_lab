load("sroots.mat");

%% initial params

x0 = [0; 0; 0; 0.2];

A = [0, 0, 1, 0;
    0, 0, 0, 1;
    -1.7117, 0, -0.3249, 0;
    0, 0, 0, -1.0004];
B = [0, 0;
    0, 0;
    0.0377, 0.0959;
    -0.1228, 0.1];
C = eye(4);
D = zeros(4,2);

Ts = 2;


%% poles
imag2ADP = [1.2982j];
ADP = [s1/Ts + imag2ADP(1,1)];
sPoles = [ADP(1,1), conj(ADP(1,1)), s2/Ts]

n = size(A,1); %size should be 4
T = min(Ts ./ (20 .* n), pi ./ (5 .* 1.2982)) %0.025

zpoles = exp(T * sPoles);

[phi, gamma] = c2d(A, B, T);

K = place(phi, gamma, zpoles);

[delta1, delta2] = rb_regsf(phi, gamma, K, T)
if(~(min(delta1, delta2) >= 0.5))
    disp('warning: stability margin below 0.5')
end
