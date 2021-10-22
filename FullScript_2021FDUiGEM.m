%开始调试
tspan=[0,150];
y0=[10,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
%[1.35812642922287,0,0,1.42680370865216,1.73512784823143,1.35818452390786,1.46410264664405,0,0,1.73512986033293,0,0,0,0,0,0,0,0,0,1.30640531526269,1.66988307956940,1.33738597942019,1.71343883384758,0,0,0,0,1.28495084599984,10,10,0.0141412486254197,0.0180757315724974,0.0154445069596673,1.06364155826607];
[T,Y]=ode45(@odefun1,tspan,y0);
%sigma70 1; sigma70lac 2; sigma70ara 3; sigma70bind 4; Msigma70 5
%sigmas 6; sigmasbind 7; sigmaslac 8; sigmasara 9; Msigmas 10
%sigma70igp2 11; sigmasigp5.7 12
%mt7p 13; t7p 14; t7plac 15
%inlac 16; outlac 29
%inara 17; outara 30
%gfp 18; mgfp 18
%p1 20; mp1 21
%p2 22; mp2 23
%ipg2 24; migp2 25
%igp57 26; migp57 27
%sigmasbind2 or sigma70bind2 28
%P3 31; MP3 32; sigma70bind3 33
%condition check 34
%----------------------------------------------------------------------------------------
GFP1=[0.930,1.243,1.511,1.920,3.744,3.819,3.826,3.357,3.461,3.938,3.491,3.513,3.528,3.238,3.260,3.446,3.72,3.61,3.7];
tGFP=12:6:120;

GFP2=[0.952917505,1.310261569,1.53360161,1.98028169,3.759557344,3.685110664,3.856338028,3.431991952,3.491549296,3.506438632,3.469215292,3.431991952,3.5138833,3.357545272,3.3277666,3.409657948,3.47665996,3.789336016,3.484104628];
GFP3=[0.752,0.945,0.990,1.280,1.913,2.583,3.559,4.206,4.526,4.288,4.593,4.497,4.482,4.645,4.288,4.534,4.564,4.571,4.340];
GFP4=[0.797,0.998,1.161,1.414,2.010,2.725,3.901,4.348,4.437,4.645,4.645,4.578,4.698,4.638,4.705,4.608,4.348,4.549,4.616];
GFP5=[0.916,1.154,1.586,2.107,2.978,4.035,4.377,4.206,4.281,4.407,4.504,4.616,4.623,4.586,4.400,4.601,4.556,4.489,4.221];
GFP6=[0.953,1.169,1.682,2.032,2.717,4.199,4.511,4.229,4.191,4.370,4.266,4.243,4.273,4.407,4.355,4.467,4.348,4.415,4.571];

subplot(4,1,1)
plot(T,Y(:,1),'.',T,Y(:,4),'-',T,Y(:,5),'-.');
legend('sigma70','sigma70bind','Msigma70')
xlabel('Time / min')
ylabel('Concentration/ mmol/ml')

subplot(4,1,2)
plot(T,Y(:,6),'-.',T,Y(:,7),'.',T,Y(:,10),'-');
legend('sigmas','sigmasbind','Msigmas')
xlabel('Time / min')
ylabel('Concentration/ mmol/ml')

subplot(4,1,3)
plot(T,Y(:,20)+Y(:,22),'.');
legend('Endogenous protein')
xlabel('Time / min')
ylabel('Concentration/ mmol/ml')

subplot(4,1,4)
plot(T,Y(:,26),'.',T,Y(:,27),T,Y(:,12));
legend('Igp5.7','MIgp5.7','sigmasIgp5.7')
xlabel('Time / min')
ylabel('Concentration/ mmol/ml')


%-------------------------------------------------------

%Igp5.7 is expressed by sigmas promoter
function dy=odefun1(t,y)

[lamdaMT7P,lamdaMGFP,lamdaMIgp2, lamdaMIgp57,lamdaMsigma70,lamdaMsigmas,lamdaMP1,lamdaMP2,lamdaMPtest]=deal(0.47/4);%unit s^-1,reference https://2019.igem.org/Team:Fudan-TSI/Model
[lamdaT7P,lamdaGFP,lamdaIgp2,lamdaIgp57,lamdasigma70,lamdasigmas,lamdaP1,lamdaP2,lamdaPtest]=deal(0.237/4);%unit s^-1 reference https://2019.igem.org/Team:Fudan-TSI/Model

%Translation Kinetic Constant
[KT7P,KGFP,KIgp2,KIgp57,Ksigmas,Ksigma70,KP1,KP2,KPtest]=deal(0.1849/4);%unit s^-1 
   
%Transcriptin Kinetic Constant of Sigma70 or Sigmas 
K70MT7P=0.585/4;  %reference https://2019.igem.org/Team:Fudan-TSI/Model [1]
K70MIgp2=0.585/4;  
KMGFP=0.5/4;  
KsMT7P=0.6/4;
KsMIgp2=0.6/4;
K70MP1=0.5847/4;
K70MPtest=0.5/4;
KsMP2=0.6/4;
K70Msigma70=0.585/4; 
K70Msigmas=0.585/4; 
KsMsigma70=0.6/4; 
KsMsigmas=0.6/4;
KsMIgp57=0.52/4;


%Release Kinetic Constant(l for lac, a for ara)
[Krlsigma70,Krasigma70, Krsigma70bind, KrlT7P,Krlsigmas, Krasigmas,Krsigmasbind, Krsigmasbind2,Krsigma70bind3]=deal(0.55/4);


K1=12/4;
K2=12/4;

condition=2.8;

%internal protein DNA片段在基因组长度比例
n1=0.05/4;
n2=0.05/4;%(1-n1*(y(16)>0)-n2*(y(17)>0))
n3=0.05/4;


%主动运输速率
Ktranslac=0.5/4;
Ktransara=0.5/4;

dy=zeros(34,1);



dy(1)=y(5)*Ksigma70+y(2)*Krlsigma70+y(4)*Krsigma70bind-y(16)*y(1)*n1*K70MT7P-y(17)*y(1)*n2*K70MIgp2-y(1)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*K70MP1-y(1)*y(24)*K1-y(1)*lamdasigma70+y(33)*Krsigma70bind3-y(1)*n3*K70MPtest;
dy(2)=y(1)*y(16)*n1*K70MT7P-y(2)*Krlsigma70;
dy(3)=y(1)*y(17)*n2*K70MIgp2-y(3)*Krasigma70;
dy(4)=y(1)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*K70MP1-y(4)*Krsigma70bind;
dy(5)=y(1)*K70Msigma70*(1-1*((y(20)+y(22))>condition))*(y(34)<1)+y(6)*(1-0.9*((y(20)>condition)*(y(34)>1)))*KsMsigma70-y(5)*lamdaMsigma70;
dy(6)=y(10)*Ksigmas+y(8)*Krlsigmas+y(9)*Krasigmas+y(7)*Krsigmasbind+y(28)*Krsigmasbind2-y(16)*y(6)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMT7P-y(17)*y(6)*n2*KsMIgp2-y(6)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMP2-y(6)*y(26)*K2-y(6)*lamdasigmas-y(6)*KsMIgp57;
dy(7)=y(6)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMP2-y(7)*Krsigmasbind;
dy(8)=y(6)*y(16)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMT7P-y(8)*Krlsigmas;
dy(9)=y(6)*y(17)*n2*KsMIgp2-y(9)*Krasigmas;
dy(10)=0.6*y(1)*K70Msigmas*((y(20)+y(22))>condition)+y(6)*KsMsigmas-y(10)*lamdaMsigmas;
dy(11)=y(1)*y(24)*K1;
dy(12)=y(6)*y(26)*K2;
dy(13)=y(1)*y(16)*n1*K70MT7P+y(6)*y(16)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMT7P-y(13)*lamdaMT7P;
dy(14)=y(13)*KT7P+y(15)*KrlT7P-y(16)*y(14)*KMGFP-y(14)*lamdaT7P;
dy(15)=y(14)*y(16)*KMGFP-y(15)*KrlT7P;
 
dy(18)=y(19)*KGFP-y(18)*lamdaGFP;
dy(19)=y(14)*y(16)*KMGFP-y(19)*lamdaMGFP;
dy(20)=y(21)*KP1-y(20)*lamdaP1;
dy(21)=y(1)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*K70MP1-y(21)*lamdaMP1;
dy(22)=y(23)*KP2-y(22)*lamdaP2;
dy(23)=y(6)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMP2-y(23)*lamdaMP2;
dy(24)=y(25)*KIgp2-y(1)*y(24)*K1-y(24)*lamdaIgp2;
dy(25)=y(1)*y(17)*n2*K70MIgp2+y(6)*y(17)*n2*KsMIgp2-y(25)*lamdaMIgp2;
dy(26)=y(27)*KIgp57-y(6)*y(26)*K2-y(26)*lamdaIgp57;
dy(27)=y(6)*y(17)*KsMIgp57-y(27)*lamdaMIgp57;
dy(28)=y(6)*KsMIgp57-y(28)*Krsigmasbind2;
 
    
dy(16)=Ktranslac*(y(29)>0)-y(16)*y(1)*n1*K70MT7P-y(16)*y(14)*KMGFP-y(6)*y(16)*(1-n1*(y(16)>0)-n2*(y(17)>0)-n3)*KsMT7P;
dy(17)=Ktransara*(y(30)>0)-y(17)*y(1)*n2*K70MIgp2-y(6)*y(17)*n2*KsMIgp2-y(6)*y(17)*KsMIgp57;
dy(29)=-Ktranslac*(y(29)>0);
dy(30)=-Ktransara*(y(30)>0);


dy(31)=y(32)*KPtest-y(31)*lamdaPtest;
dy(32)=y(1)*n3*K70MPtest-y(32)*lamdaMPtest; 
dy(33)=y(1)*n3*K70MPtest-y(33)*Krsigma70bind3;

dy(34)=1*((y(20)+y(22))>condition)*(1-(y(34)>1));

end

