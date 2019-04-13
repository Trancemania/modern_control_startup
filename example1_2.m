% run example1_1.m first
syms x
Px = poly2sym(num,x)/poly2sym(den,x);                               %symbolic transfer function with variable x
syms kpx kdx
Ts = Px*kpx/(1+Px*(kdx*x+kpx));                                     %close-loop transfer function with PD controller
invTs_bar = taylor(1/Ts, x, 0, 'Order', 4);                         %taylor series of first three terms
syms wn ksi
GMx = wn^2/(x^2+2*ksi*wn*x+wn^2);                                   %normal model transfer function
para_invGMs = coeffs(1/GMx, x, 'All');                              %symbolic coefficiets in matrix format
para_invTs_bar = coeffs(invTs_bar, x, 'All');                       %symbolic coefficiets in matrix format
[kpsol,kdsol] = solve(para_invGMs == para_invTs_bar, [kpx,kdx]);    %solve
kp = double(subs(kpsol, [wn ksi], [0.5 0.7]));                      %numeric substitution
kd = double(subs(kdsol, [wn ksi], [0.5 0.7]));                      %numeric substitution