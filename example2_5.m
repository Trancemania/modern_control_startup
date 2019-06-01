%% vc as states
syms vc(t) i(t) C x vl(t) R vr(t) L vin u
eq1 = vin == vl(t)+vr(t)+vc(t);
eq2 = vl(t) == L*diff(i(t),t);
eq3 = vr(t) == R*i(t);
eq4 = diff(vc(t),t) == diff(1/C*int(i(t),t), t);
eq4x = subs(eq4,i,x);
cancelI = solve(eq4x,x);
subi2 = subs(eq2,i,cancelI);
subi3 = subs(eq3,i,cancelI);
% rhs is available after version 2017a
% subi2 = rhs(subs(eq2,i,cancelI));
% subi3 = rhs(subs(eq3,i,cancelI));
eq = subs(subs(eq1,vl,C*L*diff(vc(t), t, t)),vr, C*R*diff(vc(t), t));
[vf,Yvar] = odeToVectorField(eq);
F_vf = matlabFunction(vf);
x_vf = F_vf(C,L,R,Yvar,u);
mat_A = equationsToMatrix(x_vf ==0, Yvar);
mat_B = equationsToMatrix(x_vf ==0, u);

%% charge as states
syms q(t)
eq4q = vc(t) == 1/C*int(i(t),t);
eq5 = diff(q(t),t) == diff(int(i(t),t),t);
eq5x = subs(eq5,i,x);
canceliq = solve(eq5x,x);
subi2q = subs(eq2,i,canceliq);
subi3q = subs(eq3,i,canceliq);
subi4q = subs(eq4q,i,canceliq);
eqq = subs(subs(subs(eq1,vl,L*diff(q(t), t, t)),vr, R*diff(q(t), t)),vc,q(t)/C);
[vfq,Yvarq] = odeToVectorField(eqq);
F_vfq = matlabFunction(vfq);
x_vfq = F_vfq(C,L,R,Yvarq,u);
mat_Aq = equationsToMatrix(x_vfq ==0, Yvarq);
mat_Bq = equationsToMatrix(x_vfq ==0, u);
fout = subs(sym('y=1/C*int(i(t),t)'),i,canceliq);

