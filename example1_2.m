syms x
Px = poly2sym(num,x)/poly2sym(den,x);
syms kpx kdx
Ts = Px*kpx/(1+Px*(kdx*x+kpx));
invTs_bar = taylor(1/Ts, x, 0, 'Order', 4);
syms wn ksi
GMx = wn^2/(x^2+2*ksi*wn*x+wn^2);
para_invGMs = coeffs(1/GMx, x, 'All');
para_invTs_bar = coeffs(invTs_bar, x, 'All');
[kpsol,kdsol] = solve(para_invGMs == para_invTs_bar, [kpx,kdx]);
kp = double(subs(kpsol, [wn ksi], [0.5 0.7]));
kd = double(subs(kdsol, [wn ksi], [0.5 0.7]));