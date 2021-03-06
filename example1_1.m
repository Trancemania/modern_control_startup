syms s y u z
eq1 = 0.5*s^2*z==u-2*(z-y)-(s*z-s*y);
eq2 = s^2*y==2*(z-y)+(s*z-s*y);         %input simultaneous equations

cancelZ = solve(eq2, z);                %cancel Z from eq2
uy = subs(eq1, z, cancelZ);             %substitute Z into eq1
Ps = subs(solve(uy, y),u,1);            %symbolic transfer function
[symnum,symden] = numden(Ps);           %extract num den
num = sym2poly(symnum);                 
den = sym2poly(symden);                 %matrix format

tf = tf(num,den);                      %transfer function
%num = cell2mat(tf.numerator(1,1));     
%den = cell2mat(tf.denominator(1,1));   %cell to matrix