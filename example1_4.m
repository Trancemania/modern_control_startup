%% MIMO state space
syms x1(t) x2(t) x3(t) x4(t) u1 u2;

eqs = [diff(x1(t),t) == x2(t),...
       0.5*diff(x2(t),t) == -2*(x1(t)-x3(t))-(x2(t)-x4(t))+u1,...
       diff(x3(t),t) == x4(t),...
       diff(x4(t),t) == 2*(x1(t)-x3(t))+(x2(t)-x4(t))+u2];
[vf,Yvar] = odeToVectorField(eqs);
rearr = matlabFunction(Yvar);
x_state_rearr = rearr(x1,x2,x3,x4);
x_state = [x1;x2;x3;x4];
[i,j] = find(x_state*transpose(1./x_state_rearr)==1);
TrM = zeros(4);
TrM(sub2ind([4,4],i,j)) = 1;

sort_vf = TrM*vf;
F_vf = matlabFunction(sort_vf);
x_vf = F_vf(Yvar,u1,u2);
mat_A = double(equationsToMatrix(x_vf ==0, TrM*Yvar));
mat_B = double(equationsToMatrix(x_vf ==0, [u1,u2]));
mat_C = [1 0 0 0; 0 0 1 0];

%% MIMO transfer function
syms y1 y2 s;
eq1 = 0.5*s^2*y1==u1-2*(y1-y2)-(s*y1-s*y2);
eq2 = s^2*y2==u2+2*(y1-y2)+(s*y1-s*y2);
cancely1 = solve(eq2, y1);
uy1 = subs(eq1, y1, cancely1);
tf_y1u1 = subs(solve(uy1, y2),[u1,u2],[1,0]);
tf_y1u2 = subs(solve(uy1, y2),[u1,u2],[0,1]);

cancely2 = solve(eq2, y2);
uy2 = subs(eq1, y2, cancely2);
tf_y2u1 = subs(solve(uy2, y1),[u1,u2],[1,0]);
tf_y2u2 = subs(solve(uy2, y1),[u1,u2],[0,1]);

mat_tf = [tf_y1u1, tf_y1u2; tf_y2u1, tf_y2u2];
