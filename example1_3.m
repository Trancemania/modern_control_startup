% run example1_1.m first
%% transfer function
syms kp1 kp2 kd1 kd2 yref
us2 = kp1*(yref-cancelZ)-kd1*s*cancelZ+kp2*(yref-y)-kd2*s*y;        %control input
close_loop2 = us2*Ps==y;                                            %close loop
Ts2 = subs(solve(close_loop2, y),yref,1);                           %transfer function
[symnum2,symden2] = numden(Ts2);                                    %extract num den
num2 = coeffs(symnum2, s, 'All');                                   
den2 = coeffs(symden2, s, 'All');                                   %coefficients of s

ss_Ps = ss(tf);

k = [-1.9, -1.71, -5.84, -4.45];
h = 7.75;

A = [0 1 0 0;-4 -2 4 2;0 0 0 1;2 1 -2 -1];
B = [0;1;0;0];
C = [0,0,1,0];