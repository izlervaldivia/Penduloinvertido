clear
clc
// load the data

A=[0 1 0 0;0 -0.0014 0.1271 0;0 0 0 1;0 -0.0025 19.1713 0]; 
B=[0;1.777;0;3.4296];
C=[1 0 0 0;0 0 1 0 ];
D=[0;0];
//trick to tackle when "D12 is not full rank"
//inputs u and e; outputs dy1 and dy2
sys = syslin("c",A, B, C);

Ap=sys.A;
Bp=sys.B;
Cp=sys.C;
Dp=sys.D;
////Evluamos observabilidad y controlabilidad
Cc = cont_mat(Ap,Bp)
rankCc=rank(Cc)
//observabilidad
O = obsv_mat(Ap, Cp)
rankO=rank(O)
/////Plotear valores singulares//////
///Valores matrices Scilab///
tr = trzeros(sys)
w = logspace(-3,3);
sv = svplot(sys,w);
//ploteo valores singulares//
scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Valores singulares planta incial","Frequency (rad/s)", "Amplitude (dB)");
//CONTROL LQR PLANTA SIN INTEGRADOR
Q_xx=diag([800 0 600 0]); //Weights on states
R_uu   = 1; //Weight on input
G1=lqr(sys,Q_xx,R_uu);

//Augment Plant with Integrators at Plant Input
[ns,nc]=size(Bp);
//ns= number of inputs; nc=number of controls
Ai=[Ap             Bp;
    0*ones(nc,ns) 0*ones(nc,nc)];

Bi=[0*ones(ns,nc); eye(nc)];
    
Ci=[Cp 0*ones(2,1)];

Di=0*ones(2,nc);

//CONTROL LQR PLANTA CON INTEGRADOR
//We use the ricatti equation for calculate de gain of the lqr controller
//for this we have  A'*X+X*A-X*B*X+C=0 for function X=riccati(A,B,C,'c','eigen')
C=diag([8000 0 2000 0 0])
rho=1;       //Cheap control recovery parameter 
                //The smaller the parameter, the better the recovery.
R = rho*eye(nc);//Control Weigthing Matrix


//now we calculate B
B=Bi*inv(R)*Bi';
A=Ai;
X=riccati(A,B,C,'c','eigen');
//the value of the gain G
G=inv(R)*Bi'*X; //GANANCIA LQR

//////////////////////////////////

sysi=syslin('c',Ai,Bi,Ci,Di);

////Evluamos observabilidad y controlabilidad planta aumentada
Cc1 = cont_mat(Ai,Bi)
rankCc1=rank(Cc1)
//observabilidad
O1 = obsv_mat(Ai, Ci)
rankO1=rank(O1)
/////Plotear valores singulares//////
///Valores matrices Scilab///
tr = trzeros(sysi)
w = logspace(-3,3);
sv = svplot(sysi,w);
//ploteo valores singulares//
scf(2);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Valores singulares planta aumentada","Frequency (rad/s)", "Amplitude (dB)");
//Design Kalman Filter//
Vd =0.1* eye (5,5)
Vn =0.1* eye (2 ,2)
H= lqe(sysi, Vd , Vn )
H=-H
I=eye(2,2);






