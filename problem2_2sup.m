[A,B]=jordan(mat_A);
[Ax,Bx]=jordan(mat_Ax);
J=A\mat_A*A;
Jx=Ax\mat_Ax*Ax;
[V,D]=eig(mat_A);
[Vx,Dx]=eig(mat_Ax);

P=A*inv(Ax);
T=[1,0,0,0;0,1,0,0;1,0,-1,0;0,1,0,-1];
N=T*inv(P);
M=P*inv(T);

P*mat_A*inv(P) == ??
inv(P)*mat_A*P == mat_Ax
T*mat_A*inv(T) == mat_Ax
inv(T)*mat_A*T == mat_Ax
N*mat_A*inv(N) == mat_A
inv(N)*mat_A*N == mat_A


P*mat_Ax*inv(P) == mat_A
inv(P)*mat_Ax*P == ??
T*mat_Ax*inv(T) == mat_A
inv(T)*mat_Ax*T == mat_A
N*mat_Ax*inv(N) == ??
inv(N)*mat_Ax*N == ??