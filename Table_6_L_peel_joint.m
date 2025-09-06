
% Title: Analytical solutions for predicting the interface tractions of micro/nano adhesive bonding joints
% Authors: Huaguo Wang, Hansong Ma*, Jingru Song, Hengxu Song, Xiaoming Liu*
% Corrsponding authors: mahs@lnm.imech.ac.cn (H. Ma), xiaomingliu@imech.ac.cn (X. Liu)

% This is a program for solving the systems of linear equations, 
% as shown in Table 6 on pages 25-26 for L-peel joint.

% Please enter the material properties, geometries, and boundary conditions.

%% 1. Balanced case - with slender adherends described by E-B theory
% characteristic roots

o = hA/l;
A = o/sqrt(2);
B = o/(2*sqrt(1+vA));
a3_EB = -A*cosh(A/2)/(2*sinh(A/2)-A*cosh(A/2));
a7_EB = -B*cosh(B/2)/(2*sinh(B/2)-B*cosh(B/2));

t1EB_Bcoeffs = [1,0,-(4* a7_EB* barGA * (1+ (1/(barE2*barh2)))/barhA),0];
t1EB_Brt = roots(t1EB_Bcoeffs);
gamma3 = t1EB_Brt(1,:); gamma1 = t1EB_Brt(2,:); gamma2 = t1EB_Brt(3,:);
t3EB_Bcoeffs = [1,0,0,0,(12* a3_EB* barEA * (1+ (1/barh2))/barhA)];
t3EB_Brt = roots(t3EB_Bcoeffs);
gamma4 = t3EB_Brt(1,:); gamma5 = t3EB_Brt(2,:);gamma6 = t3EB_Brt(3,:);gamma7 = t3EB_Brt(4,:);       

% F8 and F11 satisfy:
A1 = [1,0,-1,0
      1,6,1/(barE2*barh2),6
      0,1,0,-1/barh2
      -(1+barh2),2,0,2];
B1 = [0
      0
      0
      barM0 + barP0*barr0 - barQ1*barl1 + barQ2*barl2 ];
F8F10F11F13 = A1\B1; 

F8 = F8F10F11F13(1,:);
F11 = F8F10F11F13(3,:);

% F1 - F3 satisfy:
F3 = 3*barQ2/(2*(1+(1/barE2*barh2)));

A2 = [exp(-gamma1*bara)/gamma1, exp(-gamma2*bara)/gamma2
      exp(gamma1*bara)/gamma1, exp(gamma2*bara)/gamma2];
B2 = [F3*bara - F11
    -F3*bara - F8];
F1F2 = A2\B2 ;

F1 = F1F2(1,:);
F2 = F1F2(2,:);

% F4 - F7 satisfy:
A3 = [exp(-gamma4*bara)/gamma4, exp(-gamma5*bara)/gamma5, exp(-gamma6*bara)/gamma6, exp(-gamma7*bara)/gamma7
    exp(gamma4*bara)/gamma4, exp(gamma5*bara)/gamma5, exp(gamma6*bara)/gamma6, exp(gamma7*bara)/gamma7
    exp(-gamma4*bara)/gamma4^2, exp(-gamma5*bara)/gamma5^2, exp(-gamma6*bara)/gamma6^2, exp(-gamma7*bara)/gamma7^2
    exp(gamma4*bara)/gamma4^2, exp(gamma5*bara)/gamma5^2, exp(gamma6*bara)/gamma6^2, exp(gamma7*bara)/gamma7^2];
B3 = [-barP0 - barQ2/(1+barh2)
    -barQ2/(1+barh2)
    -(barM0*barh2 + barP0*barr0*barh2 + barQ1*barl1)/(1+barh2)
    -barQ2*barl2/(1+barh2)];
F4F5F6F7 = A3\B3; 

F4 = F4F5F6F7(1,:);
F5 = F4F5F6F7(2,:);
F6 = F4F5F6F7(3,:);
F7 = F4F5F6F7(4,:);
%% 2. Unbalanced case - with slender adherends described by E-B theory
% characteristic roots

k1 = -4*a7_EB*barGA*(1+(1/barE2*barh2))/barhA;
k2 = 12*a3_EB*barEA*(1+(1/barE2*barh2^3))/barhA;
k4 = 1/(6*a3_EB*barEA*(1-(1/barE2*barh2^2))/barhA);
t = 6*a7_EB*barGA*((1/barE2*barh2^2)-1)/barhA;
k3 = k1*k2-t/k4;

p = (3*k2 - k1^2)/3;
q = (2*k1^3 - 9*k1*k2 + 27*k3)/27;
alpha = (-q/2 + (q^2/4 + p^3/27)^0.5 )^(1/3);
beta = (-q/2 - (q^2/4 + p^3/27)^0.5 )^(1/3);
zita1 = alpha + beta;
omiga = -0.5 + (3)^0.5 /2 *1i; 
zita2 = omiga*alpha + omiga^2 * beta;
zita3 = omiga^2 *alpha + omiga*beta;

kesi1 = zita1 - k1/3;
kesi2 = zita2 - k1/3;
kesi3 = zita3 - k1/3;
lambda1 = (kesi1)^0.5; lambda2 = -(kesi1)^0.5;
lambda3 = (kesi2)^0.5; lambda4 = -(kesi2)^0.5;
lambda5 = (kesi3)^0.5; lambda6 = -(kesi3)^0.5;

% C7, C9, C12 satisfy:
A4 = [2/3*(1+(1/(barE2 *barh2))),-1,1/(barE2*barh2^2)
      0.5*(1-(1/(barE2*barh2^2))),-1,-1/(barE2*barh2^3)
      0, 1,-1];
B4 = [0
      0
      barQ2];
C7C9C12 = A4\B4;
C7 =  C7C9C12(1,:);
C9 = C7C9C12(2,:);
C12 = C7C9C12(3,:);

% C8, C10, C11, C13 satisfy:
A5 = [1,0,-1,0
      1,6,1/(barE2*barh2),6/(barE2*barh2^2)
      0,1,0,-1/(barE2*barh2^3)
      -(1+barh2),2,0,2];
B5 = [0
      0
      0
      barM0 + barP0*barr0 - barQ1*barl1 + barQ2*barl2 ];
C8C10C11C13 = A5\B5; 

C8 = C8C10C11C13(1,:);
C10 = C8C10C11C13(2,:);
C11 = C8C10C11C13(3,:);
C13 = C8C10C11C13(4,:);

% C1-C6 satisfy:
A6 = [exp(-lambda1*bara)/lambda1,exp(-lambda2*bara)/lambda2,exp(-lambda3*bara)/lambda3,exp(-lambda4*bara)/lambda4,exp(-lambda5*bara)/lambda5,exp(-lambda6*bara)/lambda6
     exp(lambda1*bara)/lambda1,exp(lambda2*bara)/lambda2,exp(lambda3*bara)/lambda3,exp(lambda4*bara)/lambda4,exp(lambda5*bara)/lambda5,exp(lambda6*bara)/lambda6
     k4*exp(-lambda1*bara)*(lambda1^2+k2/(lambda1^2)),k4*exp(-lambda2*bara)*(lambda2^2+k2/(lambda2^2)),k4*exp(-lambda3*bara)*(lambda3^2+k2/(lambda3^2)),k4*exp(-lambda4*bara)*(lambda4^2+k2/(lambda4^2)),k4*exp(-lambda5*bara)*(lambda5^2+k2/(lambda5^2)),k4*exp(-lambda6*bara)*(lambda6^2+k2/(lambda6^2))
     k4*exp(lambda1*bara)*(lambda1^2+k2/(lambda1^2)) ,k4*exp(lambda2*bara)*(lambda2^2+k2/(lambda2^2)), k4*exp(lambda3*bara)*(lambda3^2+k2/(lambda3^2)), k4*exp(lambda4*bara)*(lambda4^2+k2/(lambda4^2)) ,k4*exp(lambda5*bara)*(lambda5^2+k2/(lambda5^2)), k4*exp(lambda6*bara)*(lambda6^2+k2/(lambda6^2))
     exp(-lambda1*bara)/lambda1^2, exp(-lambda2*bara)/lambda2^2, exp(-lambda3*bara)/lambda3^2 ,exp(-lambda4*bara)/lambda4^2, exp(-lambda5*bara)/lambda5^2, exp(-lambda6*bara)/lambda6^2
     exp(lambda1*bara)/lambda1^2,  exp(lambda2*bara)/lambda2^2,  exp(lambda3*bara)/lambda3^2  ,exp(lambda4*bara)/lambda4^2,  exp(lambda5*bara)/lambda5^2,  exp(lambda6*bara)/lambda6^2];
B6 = [-barP0 - C9
      -C9
      C7*bara - C11
      -C7*bara - C8
      -(barM0*barh2 + barP0*barr0*barh2 + barQ1*barl1 - (C9*barh2 + C12)*bara - C10*barh2 +C13)/(1+barh2)
       (barQ2*barl2 - (C9*barh2 + C12)*bara + C10*barh2 -C13)/(1+barh2)];

C1C2C3C4C5C6 = A6\B6;
C1 = C1C2C3C4C5C6(1,:);C2 = C1C2C3C4C5C6(2,:);C3 = C1C2C3C4C5C6(3,:);
C4 = C1C2C3C4C5C6(4,:);C5 = C1C2C3C4C5C6(5,:);C6 = C1C2C3C4C5C6(6,:);

%% 3. Balanced case - with thick adherends described by T-E theory
% characteristic roots

o = hA/l;
a3_TE = a3_EB;
a7_TE = a7_EB;

t1TE_Bcoeffs = [1,0,-(4* a7_TE* barGA * (1+ (1/(barE2*barh2)))/barhA),0];
t1TE_Broots  = roots(t1TE_Bcoeffs);
Gamma3 = t1TE_Broots(1,:); Gamma1 = t1TE_Broots(2,:); Gamma2 = t1TE_Broots(3,:);
t3TE_Bcoeffs = [1,0,-( a3_TE* barEA* ( 2*(1+v1)+(2*(1+v2)/(barE2*barh2)) ) / (ks*barhA) ),0,(12* a3_TE* barEA * (1+ (1/barh2))/barhA)];
t3TE_Broots = roots(t3TE_Bcoeffs);
Gamma4 = t3TE_Broots(1,:); Gamma5 = t3TE_Broots(2,:); Gamma6 = t3TE_Broots(3,:); Gamma7 = t3TE_Broots(4,:);

% L8, L11 satisfy:
A7 = [1,0,-1,0
      1,6,1/(barE2*barh2),6
      0,1,0,-1/barh2
      -(1+barh2),2,0,2];
B7 = [0
      0
      0
      barM0 + barP0*barr0 - barQ1*barl1 + barQ2*barl2 ];
L8L10L11L13 = A7\B7; 

L8 = L8L10L11L13(1,:);
L11 = L8L10L11L13(3,:);

% L1 - L3 satisfy:
L3 = 3*barQ2/(2*(1+(1/barE2*barh2)));

A8 = [exp(-Gamma1*bara)/Gamma1, exp(-Gamma2*bara)/Gamma2
      exp(Gamma1*bara)/Gamma1, exp(Gamma2*bara)/Gamma2];
B8 = [L3*bara - L11
    -L3*bara - L8];
L1L2 = A8\B8 ;

L1 = L1L2(1,:);
L2 = L1L2(2,:);

% L4 - L7 satisfy:
A9 = [exp(-Gamma4*bara)/Gamma4, exp(-Gamma5*bara)/Gamma5, exp(-Gamma6*bara)/Gamma6, exp(-Gamma7*bara)/Gamma7
    exp(Gamma4*bara)/Gamma4, exp(Gamma5*bara)/Gamma5, exp(Gamma6*bara)/Gamma6, exp(Gamma7*bara)/Gamma7
    exp(-Gamma4*bara)/Gamma4^2, exp(-Gamma5*bara)/Gamma5^2, exp(-Gamma6*bara)/Gamma6^2, exp(-Gamma7*bara)/Gamma7^2
    exp(Gamma4*bara)/Gamma4^2, exp(Gamma5*bara)/Gamma5^2, exp(Gamma6*bara)/Gamma6^2, exp(Gamma7*bara)/Gamma7^2];
B9 = [-barP0 - barQ2/(1+barh2)
    -barQ2/(1+barh2)
    -(barM0*barh2 + barP0*barr0*barh2 + barQ1*barl1)/(1+barh2)
    -barQ2*barl2/(1+barh2)];

L4567 = A9\B9; 

L4 = L4567(1,:);L5 = L4567(2,:); L6 = L4567(3,:); L7 = L4567(4,:);

%% 4. Unbalanced case - with thick adherends described by T-E theory

d1 = -4*a7_TE*barGA*(1+(1/barE2*barh2))/barhA;
d2 = 6*a7_TE*barGA*((1/barE2*barh2^2)-1)/barhA;
d3 = -a3_TE*barEA*( 2*(1+v1)+(2*(1+v2)/(barE2*barh2)) ) / (ks*barhA);
d4 = 12*a3_TE*barEA*(1+(1/barE2*barh2^3))/barhA;
d5 = 6*a3_TE*barEA*(1-(1/barE2*barh2^2))/barhA;

P = (3*(d4+d1*d3) - (d1+d3)^2)/3;
Q = (2*(d1+d3)^3 - 9*(d1+d3)*(d4+d1*d3) + 27*(d1*d4-d2))/27;
Alpha = (-Q/2 + (Q^2/4 + P^3/27)^0.5 )^(1/3);
Beta = (-Q/2 - (Q^2/4 + P^3/27)^0.5 )^(1/3);
Zita1 = Alpha + Beta;
Omiga = -0.5 + (3)^0.5 /2 *1i; 
Zita2 = Omiga * Alpha + Omiga^2 * Beta;
Zita3 = Omiga^2 * Alpha + Omiga * Beta;

Kesi1 = Zita1 - (d1+d3)/3;
Kesi2 = Zita2 - (d1+d3)/3;
Kesi3 = Zita3 - (d1+d3)/3;
Lambda1 = (Kesi1)^0.5; Lambda2 = -(Kesi1)^0.5;
Lambda3 = (Kesi2)^0.5; Lambda4 = -(Kesi2)^0.5;
Lambda5 = (Kesi3)^0.5; Lambda6 = -(Kesi3)^0.5;

% H7, H9, H12 satisfy:
A10 = [2/3*(1+(1/(barE2 *barh2))),-1,1/(barE2*barh2^2)
      0.5*(1-(1/(barE2*barh2^2))),-1,-1/(barE2*barh2^3)
      0, 1,-1];
B10 = [0
      0
      barQ2];
H7H9H12 = A10\B10;
H7 = H7H9H12(1,:);
H9 = H7H9H12(2,:);
H12 = H7H9H12(3,:);

% H8, H10, H11, H13 satisfy:
A11 = [1,0,-1,0
      1,6,1/(barE2*barh2),6/(barE2*barh2^2)
      0,1,0,-1/(barE2*barh2^3)
      -(1+barh2),2,0,2];
B11 = [0
       0
       0
       barM0 + barP0*barr0 - barQ1*barl1 + barQ2*barl2 ];
H8H10H11H13 = A11\B11; 

H8 = H8H10H11H13(1,:);
H10 = H8H10H11H13(2,:);
H11 = H8H10H11H13(3,:);
H13 = H8H10H11H13(4,:);

% H1-H6 satisfy:
A12 = [exp(-Lambda1*bara)/Lambda1,exp(-Lambda2*bara)/Lambda2,exp(-Lambda3*bara)/Lambda3,exp(-Lambda4*bara)/Lambda4,exp(-Lambda5*bara)/Lambda5,exp(-Lambda6*bara)/Lambda6
     exp(Lambda1*bara)/Lambda1,exp(Lambda2*bara)/Lambda2,exp(Lambda3*bara)/Lambda3,exp(Lambda4*bara)/Lambda4,exp(Lambda5*bara)/Lambda5,exp(Lambda6*bara)/Lambda6
     (1/d5)*exp(-Lambda1*bara)*(Lambda1^2+d3+d4/(Lambda1^2)),(1/d5)*exp(-Lambda2*bara)*(Lambda2^2+d3+d4/(Lambda2^2)),(1/d5)*exp(-Lambda3*bara)*(Lambda3^2+d3+d4/(Lambda3^2)),(1/d5)*exp(-Lambda4*bara)*(Lambda4^2+d3+d4/(Lambda4^2)),(1/d5)*exp(-Lambda5*bara)*(Lambda5^2+d3+d4/(Lambda5^2)),(1/d5)*exp(-Lambda6*bara)*(Lambda6^2+d3+d4/(Lambda6^2))
     (1/d5)*exp(Lambda1*bara)*(Lambda1^2+d3+d4/(Lambda1^2)) ,(1/d5)*exp(Lambda2*bara)*(Lambda2^2+d3+d4/(Lambda2^2)), (1/d5)*exp(Lambda3*bara)*(Lambda3^2+d3+d4/(Lambda3^2)), (1/d5)*exp(Lambda4*bara)*(Lambda4^2+d3+d4/(Lambda4^2)) ,(1/d5)*exp(Lambda5*bara)*(Lambda5^2+d3+d4/(Lambda5^2)), (1/d5)*exp(Lambda6*bara)*(Lambda6^2+d3+d4/(Lambda6^2))
     exp(-Lambda1*bara)/Lambda1^2, exp(-Lambda2*bara)/Lambda2^2, exp(-Lambda3*bara)/Lambda3^2 ,exp(-Lambda4*bara)/Lambda4^2, exp(-Lambda5*bara)/Lambda5^2, exp(-Lambda6*bara)/Lambda6^2
     exp(Lambda1*bara)/Lambda1^2,  exp(Lambda2*bara)/Lambda2^2,  exp(Lambda3*bara)/Lambda3^2  ,exp(Lambda4*bara)/Lambda4^2,  exp(Lambda5*bara)/Lambda5^2,  exp(Lambda6*bara)/Lambda6^2];
B12 = [-barP0 - H9
      -H9
      H7*bara - H11
      -H7*bara - H8
      -(barM0*barh2 + barP0*barr0*barh2 + barQ1*barl1 - (H9*barh2 + H12)*bara - H10*barh2 +H13)/(1+barh2)
       (barQ2*barl2 - (H9*barh2 + H12)*bara + H10*barh2 -H13)/(1+barh2)];

H1H2H3H4H5H6 = A12\B12;
H1 = H1H2H3H4H5H6(1,:);H2 = H1H2H3H4H5H6(2,:);H3 = H1H2H3H4H5H6(3,:);
H4 = H1H2H3H4H5H6(4,:);H5 = H1H2H3H4H5H6(5,:);H6 = H1H2H3H4H5H6(6,:);
