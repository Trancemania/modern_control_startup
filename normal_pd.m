% run example1_1.m first
syms kpp kdd
eq3 = y == Ps*kpp*(yref-y)+kdd*s*(yref-y);
tf_pd = subs(solve(eq3, y),yref,1);
[sym_pd_num,sym_pd_den] = numden(tf_pd);
pd_num = coeffs(sym_pd_num, s, 'All');
pd_den = coeffs(sym_pd_den, s, 'All');

tl4 = taylor(1/tf_pd, s, 0, 'Order', 4);
tl5 = taylor(1/tf_pd, s, 0, 'Order', 5);
tl6 = taylor(1/tf_pd, s, 0, 'Order', 6);

