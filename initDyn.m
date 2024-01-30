%% This m-file initiats MatSim1.0
clc
clear all
close all

%% Define Variables
global_var

%% User Must Set This pass path based on his PC
% addpath('F:\PHD_works\Joint works\Amit Kumar\MatSim project\Matpower4.1')
addpath('C:Users\EBUKA\Desktop\matpower4.1\MatSim1.0\Three_Machine\AEO')


pf = runpf(case9_1);

Vm=pf.bus(:,8);
Va=pf.bus(:,9);
Vb=Vm.*exp(1j*(Va*pi/180));

sys_data = case9_1;

baseMVA=sys_data.baseMVA;
bus=sys_data.bus;
branch1=sys_data.branch;
BRX=branch1(:,1:5);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch1);
YB=full(Ybus);
clear Ybus Yf Yt

v=Vb;

Pg=pf.gen(:,2)/baseMVA;
Qg=pf.gen(:,3)/baseMVA;
Pd=sys_data.bus(:,3)/baseMVA;
Qd=sys_data.bus(:,4)/baseMVA;
Ybus=YB;

%% Enter Dynamic parameters for the Generators
m=3;
n=9;

%         H      Xd     Xpd      Xq     Xpq    Tpd0    Tpq0
GenData=[ 23.63 0.146  0.0608  0.0969  0.0969  8.96     0.31
    6.4   0.8958 0.1198  0.8645  0.1969  6      0.535
    3.01  1.3125 0.1813  1.2578  0.25  5.89     0.6];

%        KA  TA   KE   TE     KF       TF
ExData=[ 20  0.2  1    0.314  0.063    0.35
    20  0.2  1    0.314  0.063    0.35
    20  0.2  1    0.314  0.063    0.35];

%%
H=GenData(:,1);
Xd=GenData(:,2);
Xpd=GenData(:,3);
Xq=GenData(:,4);
Xpq=GenData(:,5);
Tpd0=GenData(:,6);
Tpq0=GenData(:,7);
ws=377;
M=2*H/ws;
D(1,1)=0.1*M(1);
D(2,1)=0.2*M(2);
D(3,1)=0.3*M(3);
D=0*D;
KA=ExData(:,1);
TA=ExData(:,2);
KE=ExData(:,3);
TE=ExData(:,4);
KF=ExData(:,5);
TF=ExData(:,6);
Rs=[0 0 0];

YN=Ybus;

%% apply fault on the system
% fault=1 ; the system with fault
% fault=0 ; the system without fault

fault=1;

%% Select fault Type
% ftype=1; small change in references
% ftype=2; symetrical three phase faults
% ftype=3; load changes

ftype=2;


fPd=Pd;
fQd=Qd;

if ftype == 3
    fPd(5,1)=Pd(5,1)/2;
    fQd(5,1)=Qd(5,1)/2;
end
event=[0 1 1.10 10;
    0 1 0   0];
small_d=[
    0  1 1.1 10
    1  1.00 1 1;];
fYN = YN;
for i=1:n
    YL(i,i) = (Pd(i)-sqrt(-1)*Qd(i))/abs(v(i))^2;
    fYL(i,i) = (fPd(i)-sqrt(-1)*fQd(i))/abs(v(i))^2;
    
    YN(i,i) = YB(i,i)+YL(i,i);
    fYN(i,i) = YB(i,i)+fYL(i,i);
end
YGG = YN(1:m,1:m);
YGL = YN(1:m,m+1:n);
YLG = YN(m+1:n,1:m);
YLL = YN(m+1:n,m+1:n);

Yred = YGG-YGL*inv(YLL)*YLG;
Zred = inv(Yred);
iYLL = inv(YLL);
if fault 
    if ftype == 1   % small change in references
        fYred=Yred;
        fZred = inv(fYred);
        small_d=[
            0  1 1.1 10
            1  1.05 1 1;];
    end
    if ftype == 2 % symetrical three phase faults
        fYred=YN(1:m,1:m)-YN(1:m,m+1:n-1)*inv(YN(m+1:n-1,m+1:n-1))*YN(m+1:n-1,1:m);
        fZred = inv(fYred);
    end
    if ftype == 3 % load changes
        fYGG = fYN(1:m,1:m);
        fYGL = fYN(1:m,m+1:n);
        fYLG = fYN(m+1:n,1:m);
        fYLL = fYN(m+1:n,m+1:n);
        
        fYred = fYGG-fYGL*inv(fYLL)*fYLG;
        fZred = inv(fYred);
        ifYLL = inv(fYLL);
    end
    if ftype == 4 % line outage
        
    end
end

YLV=-inv(YN(m+1:n,m+1:n))*YN(m+1:n,1:m);
%Initial Condition Calculation----------------------------------------
V0=abs(v);
Teta0=angle(v);

for i=1:3
    IG0(i)=(Pg(i)-sqrt(-1)*Qg(i))/conj(v(i));
    delta0(i)=angle(v(i)+(Rs(i)+sqrt(-1)*Xq(i))*IG0(i));
    
    Id0(i)=real(abs(IG0(i))*exp(sqrt(-1)*(pi/2+angle(IG0(i))-delta0(i))));
    Iq0(i)=imag(abs(IG0(i))*exp(sqrt(-1)*(pi/2+angle(IG0(i))-delta0(i))));
    Vd0(i)=real(V0(i)*exp(sqrt(-1)*(pi/2+Teta0(i)-delta0(i))));
    Vq0(i)=imag(V0(i)*exp(sqrt(-1)*(pi/2+Teta0(i)-delta0(i))));
    
    Epq0(i)=Vq0(i)+Rs(i)*Iq0(i)+Xpd(i)*Id0(i);
    Epd0(i)=(Xq(i)-Xpq(i))*Iq0(i);
    Efd0(i)=Epq0(i)+(Xd(i)-Xpd(i))*Id0(i);
    TM(i)=Epq0(i)*Iq0(i)+(Xq(i)-Xpd(i))*Id0(i)*Iq0(i);
    VR0(i)=(KE(i)+0.0039*exp(1.555*Efd0(i)))*Efd0(i);
    Rf0(i)=KF(i)*Efd0(i)/TF(i);
    Vref(i)=V0(i)+VR0(i)/KA(i);
    Vref1(i)=V0(i)+Efd0(i)/KA(i);
end
%%   ALgebraic equation

Rsm=diag(Rs);
Xpqm=diag(Xpq);
Xpdm=diag(Xpd);

XXpdq=[Xpq-Xpd];
XXdpd=[Xd-Xpd];
XXqpq=[Xq-Xpq];

Zdq=[Rsm -Xpqm;Xpdm Rsm];
iZdq=inv(Zdq);

MT = sum(M);
