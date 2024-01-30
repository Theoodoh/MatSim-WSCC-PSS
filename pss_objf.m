function J = pss_objf(x,FunIndex,Dim)

load('sys_IO','f11')
As = f11.a;
Bs = f11.b;
Cs = f11.c;
Ds = f11.d;

Tw = 10;
KG2 = x(1);
T21 = x(3);
T22 = x(4);
T23 = x(5);
T24 = x(6);
Kpss2 = KG2*T21*T23/(T22*T24);

KG3 = x(2);
T31 = x(7);
T32 = x(8);
T33 = x(9);
T34 = x(10);
Kpss3 =  KG3*T31*T33/(T32*T34);


b2 = [KG2*T21*T23*Tw (KG2*T21*Tw + KG2*T23*Tw) KG2*Tw 0];
a2 = [T22*T24*Tw  (T22*T24 + T22*Tw + T24*Tw) (T22 + T24 + Tw) 1];

b3 = [KG3*T31*T33*Tw (KG3*T31*Tw + KG3*T33*Tw) KG3*Tw 0];
a3 = [T32*T34*Tw  (T32*T34 + T32*Tw + T34*Tw) (T32 + T34 + Tw) 1];

[A_1 B_1 C_1 D_1 ]= tf2ss(b2,a2);
[A_2 B_2 C_2 D_2 ]= tf2ss(b3,a3);

Af = blkdiag(A_1,A_2);
Bf = blkdiag(B_1,B_2);
Cf = blkdiag(C_1,C_2);
Df = blkdiag(D_1,D_2);

Asys_1 = As + Bs*Df*Cs;
Asys_2 = Bs*Cf;
Asys_3 = Bf*Cs;
Asys_4 = Af + Bf*Ds*Cf;
Asys = [Asys_1 Asys_2;
    Asys_3 Asys_4];

egs = eig(Asys);

[z_val z_idx]=sort(abs(egs),'descend');
egs_new=egs;
egs_new(z_idx(end-1:end))=[];

%% unstable modes
ss_idx = find(real(egs_new)>0);
uss = egs_new(ss_idx);

%% EM modes
% Damp=-real(egs)./sqrt(real(egs).^2+imag(egs).^2)
freq = abs(imag(egs_new))/(2*pi);
em_idx = find(freq>0 & freq<3);

objf = max(real(egs_new(em_idx)))+sum(egs_new(ss_idx));

if isempty(objf)
    objf = max(real(egs_new));
end
J = objf;