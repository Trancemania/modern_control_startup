syms x1(t) x2(t) x3(t) x4(t) u;
syms z1(t) z2(t);
eq1 = 0.5*diff(z1(t),t,t) == u-2*(z1(t)-z2(t))-(diff(z1(t),t)-diff(z2(t),t));
eq2 = diff(z2(t),t,t) == 2*(z1(t)-z2(t))+(diff(z1(t),t)-diff(z2(t),t));

[vf,Yvar] = odeToVectorField([eq1,eq2]);

rearr = matlabFunction(Yvar);
z_state_rearr = rearr(diff(z1,t),diff(z2,t),z1,z2);%rearr=@(Dz1,Dz2,z1,z2)...
z_state = [z1;diff(z1,t);z2;diff(z2,t)];
[i,j] = find(z_state*transpose(1./z_state_rearr)==1);
TrM = zeros(4);
TrM(sub2ind([4,4],i,j)) = 1;
sort_vf = TrM*vf;
F_vf = matlabFunction(sort_vf);
z_vf = F_vf(Yvar,u); 
mat_A = double(equationsToMatrix(z_vf ==0, TrM*Yvar));
mat_B = double(equationsToMatrix(z_vf ==0, u));

syms x11 x22 x33 x44 z11 z22 dz11 dz22
dz_state = [z11; dz11; z22; dz22];
fun = formula(mat_A*dz_state+mat_B*u);

dz1 = fun(1);
ddz1 = fun(2);
dz2 = fun(3);
ddz2 = fun(4);
dx1 = dz1;
dx2 = ddz1;
dx3 = dz1-dz2;
dx4 = ddz1-ddz2;

funxz = [x11 == z11, x22 == dz11, x33 == z11 - z22, x44 == dz11 - dz22];
[cancelz11,canceldz11,cancelz22,canceldz22]=solve(funxz,[z11,dz11,z22,dz22]);
rhsfunx = subs([dx1,dx2,dx3,dx4],[z11,dz11,z22,dz22],[cancelz11,canceldz11,cancelz22,canceldz22]);
rhsfunxt = subs(rhsfunx,[x11,x22,x33,x44],[x1,x2,x3,x4]);
odext = formula([diff(x1,t),diff(x2,t),diff(x3,t),diff(x4,t)])==rhsfunxt;
[vfx,Yvarx] = odeToVectorField(odext);

rearrx = matlabFunction(Yvarx);
x_state_rearr = rearrx(x1,x2,x3,x4);
x_state = [x1;x2;x3;x4];
[ix,jx] = find(x_state*transpose(1./x_state_rearr)==1);
TrMx = zeros(4);
TrMx(sub2ind([4,4],ix,jx)) = 1;
sort_vfx = TrMx*vfx;
F_vfx = matlabFunction(sort_vfx);
x_vf = F_vfx(Yvarx,u); 
mat_Ax = double(equationsToMatrix(x_vf ==0, TrMx*Yvarx));
mat_Bx = double(equationsToMatrix(x_vf ==0, u));


Y = sym('Y',[4 4]);
soltm = solve(Y*[z11;dz11;z22;dz22]==[z11;dz11;z11-z22;dz11-dz22],Y);
struct2cell
reshape

Y = sym('Y',[4 4]);
assume(Y,'integer');
syms a b c d;
sol = solve(Y*[a;b;c;d]==[a;b;a-c;b-d],Y);
