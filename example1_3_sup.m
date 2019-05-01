% this supplyment shows how to generate state space from differential
% equations
syms x1(t) x2(t) x3(t) x4(t) u;

eqs = [diff(x1(t),t) == x2(t),...
       0.5*diff(x2(t),t) == -2*(x1(t)-x3(t))-(x2(t)-x4(t))+u,...
       diff(x3(t),t) == x4(t),...
       diff(x4(t),t) == 2*(x1(t)-x3(t))+(x2(t)-x4(t))];
[vf,Yvar] = odeToVectorField(eqs);                                  % ODE to vector field

rearr = matlabFunction(Yvar);
x_state_rearr = rearr(x1,x2,x3,x4);
x_state = [x1;x2;x3;x4];
[i,j] = find(x_state*transpose(1./x_state_rearr)==1);
TrM = zeros(4);
TrM(sub2ind([4,4],i,j)) = 1;                                        %since odeToVectorField rearrange states, this
                                                                    %step gets transformation matrix to arrange 
                                                                    %correctly
sort_vf = TrM*vf;
F_vf = matlabFunction(sort_vf);
x_vf = F_vf(Yvar,u);
mat_A = double(equationsToMatrix(x_vf ==0, TrM*Yvar));
mat_B = double(equationsToMatrix(x_vf ==0, u));

y = x3;
[ii,jj] = find(y*transpose(1./x_state)==1);
mat_C = zeros(1,4);
mat_C(sub2ind([1,4],ii,jj)) = 1;
sys = ss(mat_A,mat_B,mat_C,0);