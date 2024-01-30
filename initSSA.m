%% this m-file perfrmes small-signal analysis
initDyn


%  AEO
x = [3.9733      2.7192     0.57378    0.020083     0.50897    0.035874     0.73904    0.026834     0.99401    0.030678];

Tw = 10;
KG2 = x(1);
T21 = x(3);
T22 = x(4);
T23 = x(5);
T24 = x(6);
Kpss2=KG2*T21*T23/(T22*T24);

KG3 = x(2);
T31 = x(7);
T32 = x(8);
T33 = x(9);
T34 = x(10);
Kpss3=KG3*T31*T33/(T32*T34);

%% Linearize Power System
% f11=linmod('sys3m');
% f11=linmod('sys3m_IO');
f11=linmod('sys3m_pss');

% dx/dt = A.x + B.u
% y = C.x + D.u
Asys=f11.a;
Bsys=f11.b;
Csys=f11.c;
Dsys=f11.d;

%% Calculate Eigenvalues
egs = eig(Asys)
Ns=length(egs);

Damp=-real(egs)./sqrt(real(egs).^2+imag(egs).^2)
freq=abs(imag(egs))/(2*pi)


%% calculae Participation Factors
[Vs,D_eig] = eig(Asys);
Ws=inv(Vs);
for i=1:Ns
    for k=1:Ns
        Pfact1(k,i)=abs(Vs(k,i))*abs(Ws(i,k));
    end
end

for i=1:Ns
     Pfact(i,:)=Pfact1(i,:)/sum(Pfact1(i,:));
end

for i=1:Ns
    [s_val s_idx]=sort(Pfact(:,i),'descend');
    mod_idx(i,:)=s_idx(1:4)';
    pf_fact(i,:)=s_val(1:4)';
end
mod_idx;
pf_fact;