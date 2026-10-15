% Dynamic Identification Model (Newton-Euler method)
% 
% Geometric parameters
% j       ant     sigma   mu      gamma   b       alpha   d       theta       r       
% 1       0       0       1       0       0       0       0       th1         0       
% 2       1       0       1       0       0       pi/2    D2      th2 + pi/2  0       
% 3       2       0       1       0       0       0       D3      th3         0       
% 4       3       0       1       0       0       pi/2    D4      th4        RL4     
% 5       4       0       1       0       0       -pi/2   0       th5         0       
% 6       5       0       1       0       0       pi/2    0       th6        RL6     
% 
% Dynamic inertia parameters
% j       XX      XY      XZ      YY      YZ      ZZ      MX      MY      MZ      M       IA      
% 1       0       0       0       0       0       ZZR1    0       0       0       0       0       
% 2       XXR2    XY2     XZR2    0       YZ2     ZZR2    MXR2    MY2     0       0       0       
% 3       XXR3    XYR3    XZ3     0       YZ3     ZZR3    MXR3    MYR3    0       0       IA3     
% 4       XXR4    XY4     XZ4     0       YZ4     ZZR4    MX4     MYR4    0       0       IA4     
% 5       XXR5    XY5     XZ5     0       YZ5     ZZR5    MX5     MYR5    0       0       IA5     
% 6       XXR6    XY6     XZ6     0       YZ6     ZZ6     MX6     MY6     0       0       IA6     
% 
% External forces and joint parameters
% j       FX      FY      FZ      CX      CY      CZ      FS      FV      QP      QDP     GAM     eta     k       
% 1       0       0       0       0       0       0       FS1     FV1     QP1     QDP1    GAM1    0       0       
% 2       0       0       0       0       0       0       FS2     FV2     QP2     QDP2    GAM2    0       0       
% 3       0       0       0       0       0       0       FS3     FV3     QP3     QDP3    GAM3    0       0       
% 4       0       0       0       0       0       0       FS4     FV4     QP4     QDP4    GAM4    0       0       
% 5       0       0       0       0       0       0       FS5     FV5     QP5     QDP5    GAM5    0       0       
% 6       0       0       0       0       0       0       FS6     FV6     QP6     QDP6    GAM6    0       0       
% 
% Base velicities parameters
% axis    W0      WP0     V0      VP0     G       
% X       0       0       0       0       0       
% Y       0       0       0       0       0       
% Z       0       0       0       0       G3      
% 
% Equations:

function w_min = min_regressor(q, dq, ddq,D2,D3,D4,RL4,RL6) %unit: rad. m 

    th1=q(1);
    th2=q(2);
    th3=q(3);
    th4=q(4);
    th5=q(5);
    th6=q(6);
    
    QP1=dq(1);
    QP2=dq(2);
    QP3=dq(3);
    QP4=dq(4);
    QP5=dq(5);
    QP6=dq(6);
    
    QDP1=ddq(1);
    QDP2=ddq(2);
    QDP3=ddq(3);
    QDP4=ddq(4);
    QDP5=ddq(5);
    QDP6=ddq(6);

    %% PARAMETERS
    % force/torque acting on ee
    % fe1 = 0; fe2 = 0; fe3 = 0; ne1 = 0; ne2 = 0; ne3 = 0;
    % modified DH parameters
    
    % gravitional acceleration
    G3 = 9.81;   % in m
    K1=153;
    K2=81;
    K3=-121;
    K4=-91.8;
    K5=-81;
    K6=-45.6875;
    % pnum_min = evalin('base', 'pnum_min');
    w_min=zeros(6,62); % number of base inertia parameters
    %% matrix
    C1 = cos(th1);
    S1 = sin(th1);
    C2 = sin(th2);
    S2 = cos(th2);
    C3 = cos(th3);
    S3 = sin(th3);
    C4 = cos(th4);
    S4 = sin(th4);
    C5 = cos(th5);
    S5 = sin(th5);
    C6 = cos(th6);
    S6 = sin(th6);
    DV61 = QP1^2;
    W12 = QP1*S2;
    W22 = -C2*QP1;
    WP12 = QDP1*S2 + QP2*W22;
    WP22 = -C2*QDP1 - QP2*W12;
    DV22 = W12*W22;
    DV32 = QP2*W12;
    DV42 = W22^2;
    DV52 = QP2*W22;
    DV62 = QP2^2;
    U112 = -DV42 - DV62;
    U212 = DV22 + QDP2;
    U312 = DV32 - WP22;
    U322 = DV52 + WP12;
    VSP12 = -D2*DV61;
    VSP22 = D2*QDP1;
    VP12 = -C2*VSP12 - G3*S2;
    VP22 = C2*G3 - S2*VSP12;
    W13 = C3*W12 + S3*W22;
    W23 = C3*W22 - S3*W12;
    W33 = QP2 + QP3;
    WP13 = C3*WP12 + QP3*W23 + S3*WP22;
    WP23 = C3*WP22 - QP3*W13 - S3*WP12;
    WP33 = QDP2 + QDP3;
    DV13 = W13^2;
    DV23 = W13*W23;
    DV33 = W13*W33;
    DV43 = W23^2;
    DV53 = W23*W33;
    DV63 = W33^2;
    U113 = -DV43 - DV63;
    U213 = DV23 + WP33;
    U313 = DV33 - WP23;
    U123 = DV23 - WP33;
    U223 = -DV13 - DV63;
    U323 = DV53 + WP13;
    VSP13 = D3*U112 + VP12;
    VSP23 = D3*U212 + VP22;
    VSP33 = D3*U312 - VSP22;
    VP13 = C3*VSP13 + S3*VSP23;
    VP23 = C3*VSP23 - S3*VSP13;
    W14 = C4*W13 + S4*W33;
    W24 = C4*W33 - S4*W13;
    W34 = QP4 - W23;
    WP14 = C4*WP13 + QP4*W24 + S4*WP33;
    WP24 = C4*WP33 - QP4*W14 - S4*WP13;
    WP34 = QDP4 - WP23;
    DV14 = W14^2;
    DV24 = W14*W24;
    DV34 = W14*W34;
    DV44 = W24^2;
    DV54 = W24*W34;
    DV64 = W34^2;
    U114 = -DV44 - DV64;
    U214 = DV24 + WP34;
    U314 = DV34 - WP24;
    U124 = DV24 - WP34;
    U224 = -DV14 - DV64;
    U324 = DV54 + WP14; 
    VSP14 = D4*U113 - RL4*U123 + VP13;
    VSP24 = D4*U213 - RL4*U223 + VP23;
    VSP34 = D4*U313 - RL4*U323 + VSP33;
    VP14 = C4*VSP14 + S4*VSP34;
    VP24 = C4*VSP34 - S4*VSP14;
    W15 = C5*W14 - S5*W34;
    W25 = -C5*W34 - S5*W14;
    W35 = QP5 + W24;
    WP15 = C5*WP14 + QP5*W25 - S5*WP34;
    WP25 = -C5*WP34 - QP5*W15 - S5*WP14;
    WP35 = QDP5 + WP24;
    DV15 = W15^2;
    DV25 = W15*W25;
    DV35 = W15*W35;
    DV45 = W25^2;
    DV55 = W25*W35;
    DV65 = W35^2;
    U115 = -DV45 - DV65;
    U215 = DV25 + WP35;
    U315 = DV35 - WP25;
    U125 = DV25 - WP35;
    U225 = -DV15 - DV65;
    U325 = DV55 + WP15;
    VP15 = C5*VP14 + S5*VSP24;
    VP25 = C5*VSP24 - S5*VP14;
    W16 = C6*W15 + S6*W35;
    W26 = C6*W35 - S6*W15;
    W36 = QP6 - W25;
    WP16 = C6*WP15 + QP6*W26 + S6*WP35;
    WP26 = C6*WP35 - QP6*W16 - S6*WP15;
    WP36 = QDP6 - WP25;
    DV16 = W16^2;
    DV26 = W16*W26;
    DV36 = W16*W36;
    DV46 = W26^2;
    DV56 = W26*W36;
    DV66 = W36^2;
    U116 = -DV46 - DV66;
    U216 = DV26 + WP36;
    U316 = DV36 - WP26;
    U126 = DV26 - WP36;
    U226 = -DV16 - DV66;
    U326 = DV56 + WP16;
    VSP16 = -RL6*U125 + VP15;
    VSP26 = -RL6*U225 + VP25;
    VSP36 = -RL6*U325 + VP24;
    VP16 = C6*VSP16 + S6*VSP36;
    VP26 = C6*VSP36 - S6*VSP16;
    DG1ZZR1 = QDP1;
    % DG1FS1 = sign(QP1);
    % DG1FV1 = QP1;
    N13XXR2 = -C2*DV32 + S2*WP12;
    DG1XXR2 = N13XXR2;
    DG2XXR2 = -DV22;
    No21XY2 = -QP2*W12 + WP22;
    No22XY2 = QP2*W22 + WP12;
    No23XY2 = W12^2 - W22^2;    
    N13XY2 = -C2*No22XY2 + No21XY2*S2;
    DG1XY2 = N13XY2;
    DG2XY2 = No23XY2;
    No21XZR2 = QDP2 + W12*W22;
    No22XZR2 = QP2^2 - W12^2;
    No23XZR2 = -QP2*W22 + WP12;   
    N13XZR2 = -C2*No22XZR2 + No21XZR2*S2;
    DG1XZR2 = N13XZR2;
    DG2XZR2 = No23XZR2;
    No21YZ2 = -QP2^2 + W22^2;
    No22YZ2 = QDP2 - W12*W22;
    No23YZ2 = QP2*W12 + WP22;  
    N13YZ2 = -C2*No22YZ2 + No21YZ2*S2;
    DG1YZ2 = N13YZ2;
    DG2YZ2 = No23YZ2;
    N13ZZR2 = C2*DV32 + DV52*S2;
    DG1ZZR2 = N13ZZR2;
    DG2ZZR2 = QDP2;
    N13MXR2 = -C2*VSP22 - D2*U312;
    DG1MXR2 = N13MXR2;
    DG2MXR2 = VP22;
    N13MY2 = -D2*U322 - S2*VSP22;
    DG1MY2 = N13MY2;
    DG2MY2 = -VP12;
    DG2FS2 = sign(QP2);
    DG2FV2 = QP2;
    N21XXR3 = C3*WP13 - DV33*S3;
    N22XXR3 = C3*DV33 + S3*WP13;    
    N13XXR3 = -C2*N22XXR3 + N21XXR3*S2;
    DG1XXR3 = N13XXR3;
    DG2XXR3 = -DV23;
    DG3XXR3 = -DV23;
    No31XYR3 = -W13*W33 + WP23;
    No32XYR3 = W23*W33 + WP13;
    No33XYR3 = W13^2 - W23^2;
    N21XYR3 = C3*No31XYR3 - No32XYR3*S3;
    N22XYR3 = C3*No32XYR3 + No31XYR3*S3;    
    N13XYR3 = -C2*N22XYR3 + N21XYR3*S2;
    DG1XYR3 = N13XYR3;
    DG2XYR3 = No33XYR3;
    DG3XYR3 = No33XYR3;
    No31XZ3 = W13*W23 + WP33;
    No32XZ3 = -W13^2 + W33^2;
    No33XZ3 = -W23*W33 + WP13;
    N21XZ3 = C3*No31XZ3 - No32XZ3*S3;
    N22XZ3 = C3*No32XZ3 + No31XZ3*S3;    
    N13XZ3 = -C2*N22XZ3 + N21XZ3*S2;
    DG1XZ3 = N13XZ3;
    DG2XZ3 = No33XZ3;
    DG3XZ3 = No33XZ3;
    No31YZ3 = W23^2 - W33^2;
    No32YZ3 = -W13*W23 + WP33;
    No33YZ3 = W13*W33 + WP23;
    N21YZ3 = C3*No31YZ3 - No32YZ3*S3;
    N22YZ3 = C3*No32YZ3 + No31YZ3*S3;   
    N13YZ3 = -C2*N22YZ3 + N21YZ3*S2;
    DG1YZ3 = N13YZ3;
    DG2YZ3 = No33YZ3;
    DG3YZ3 = No33YZ3;
    N21ZZR3 = C3*DV53 + DV33*S3;
    N22ZZR3 = -C3*DV33 + DV53*S3;    
    N13ZZR3 = -C2*N22ZZR3 + N21ZZR3*S2;
    DG1ZZR3 = N13ZZR3;
    DG2ZZR3 = WP33;
    DG3ZZR3 = WP33;    
    FDI32MXR3 = C3*U213 + S3*U113;
    N21MXR3 = S3*VSP33;
    N22MXR3 = -C3*VSP33 - D3*U313;
    N23MXR3 = D3*FDI32MXR3 + VP23;
    N13MXR3 = -C2*N22MXR3 - D2*U313 + N21MXR3*S2;
    DG1MXR3 = N13MXR3;
    DG2MXR3 = N23MXR3;
    DG3MXR3 = VP23;  
    FDI32MYR3 = C3*U223 + S3*U123;
    N21MYR3 = C3*VSP33;
    N22MYR3 = -D3*U323 + S3*VSP33;
    N23MYR3 = D3*FDI32MYR3 - VP13;
    N13MYR3 = -C2*N22MYR3 - D2*U323 + N21MYR3*S2;
    DG1MYR3 = N13MYR3;
    DG2MYR3 = N23MYR3;
    DG3MYR3 = -VP13;
    DG3IA3 = K3^2*QDP3;  %*
    DG3FS3 = sign(QP3);
    DG3FV3 = QP3;
    N31XXR4 = C4*WP14 - DV34*S4;
    N33XXR4 = C4*DV34 + S4*WP14;
    N21XXR4 = C3*N31XXR4 - DV24*S3;
    N22XXR4 = C3*DV24 + N31XXR4*S3;    
    N13XXR4 = -C2*N22XXR4 + N21XXR4*S2;
    DG1XXR4 = N13XXR4;
    DG2XXR4 = N33XXR4;
    DG3XXR4 = N33XXR4;
    DG4XXR4 = -DV24;
    No41XY4 = -W14*W34 + WP24;
    No42XY4 = W24*W34 + WP14;
    No43XY4 = W14^2 - W24^2;
    N31XY4 = C4*No41XY4 - No42XY4*S4;
    N33XY4 = C4*No42XY4 + No41XY4*S4;
    N21XY4 = C3*N31XY4 + No43XY4*S3;
    N22XY4 = -C3*No43XY4 + N31XY4*S3;   
    N13XY4 = -C2*N22XY4 + N21XY4*S2;
    DG1XY4 = N13XY4;
    DG2XY4 = N33XY4;
    DG3XY4 = N33XY4;
    DG4XY4 = No43XY4;
    No41XZ4 = W14*W24 + WP34;
    No42XZ4 = -W14^2 + W34^2;
    No43XZ4 = -W24*W34 + WP14;
    N31XZ4 = C4*No41XZ4 - No42XZ4*S4;
    N33XZ4 = C4*No42XZ4 + No41XZ4*S4;
    N21XZ4 = C3*N31XZ4 + No43XZ4*S3;
    N22XZ4 = -C3*No43XZ4 + N31XZ4*S3;   
    N13XZ4 = -C2*N22XZ4 + N21XZ4*S2;
    DG1XZ4 = N13XZ4;
    DG2XZ4 = N33XZ4;
    DG3XZ4 = N33XZ4;
    DG4XZ4 = No43XZ4;
    No41YZ4 = W24^2 - W34^2;
    No42YZ4 = -W14*W24 + WP34;
    No43YZ4 = W14*W34 + WP24;
    N31YZ4 = C4*No41YZ4 - No42YZ4*S4;
    N33YZ4 = C4*No42YZ4 + No41YZ4*S4;
    N21YZ4 = C3*N31YZ4 + No43YZ4*S3;
    N22YZ4 = -C3*No43YZ4 + N31YZ4*S3;    
    N13YZ4 = -C2*N22YZ4 + N21YZ4*S2;
    DG1YZ4 = N13YZ4;
    DG2YZ4 = N33YZ4;
    DG3YZ4 = N33YZ4;
    DG4YZ4 = No43YZ4;
    N31ZZR4 = C4*DV54 + DV34*S4;
    N33ZZR4 = -C4*DV34 + DV54*S4;
    N21ZZR4 = C3*N31ZZR4 + S3*WP34;
    N22ZZR4 = -C3*WP34 + N31ZZR4*S3; 
    N13ZZR4 = -C2*N22ZZR4 + N21ZZR4*S2;
    DG1ZZR4 = N13ZZR4;
    DG2ZZR4 = N33ZZR4;
    DG3ZZR4 = N33ZZR4;
    DG4ZZR4 = WP34;
    FDI41MX4 = C4*U114 - S4*U214;
    FDI43MX4 = C4*U214 + S4*U114;
    N31MX4 = -FDI43MX4*RL4 - S4*VSP24;
    N32MX4 = -D4*FDI43MX4 - VP24;
    N33MX4 = C4*VSP24 - D4*U314 + FDI41MX4*RL4;  
    FDI32MX4 = -C3*U314 + FDI41MX4*S3;
    N21MX4 = C3*N31MX4 - N32MX4*S3;
    N22MX4 = C3*N32MX4 - D3*FDI43MX4 + N31MX4*S3;
    N23MX4 = D3*FDI32MX4 + N33MX4;
    N13MX4 = -C2*N22MX4 - D2*FDI43MX4 + N21MX4*S2;
    DG1MX4 = N13MX4;
    DG2MX4 = N23MX4;
    DG3MX4 = N33MX4;
    DG4MX4 = VP24;
    FDI41MYR4 = C4*U124 - S4*U224;
    FDI43MYR4 = C4*U224 + S4*U124;
    N31MYR4 = -C4*VSP24 - FDI43MYR4*RL4;
    N32MYR4 = -D4*FDI43MYR4 + VP14;
    N33MYR4 = -D4*U324 + FDI41MYR4*RL4 - S4*VSP24;  
    FDI32MYR4 = -C3*U324 + FDI41MYR4*S3;
    N21MYR4 = C3*N31MYR4 - N32MYR4*S3;
    N22MYR4 = C3*N32MYR4 - D3*FDI43MYR4 + N31MYR4*S3;
    N23MYR4 = D3*FDI32MYR4 + N33MYR4;
    N13MYR4 = -C2*N22MYR4 - D2*FDI43MYR4 + N21MYR4*S2;
    DG1MYR4 = N13MYR4;
    DG2MYR4 = N23MYR4;
    DG3MYR4 = N33MYR4;
    DG4MYR4 = -VP14;
    DG4IA4 = K4^2*QDP4; %*
    % DG4FS4 = sign(QP4);
    % DG4FV4 = QP4;
    N41XXR5 = C5*WP15 - DV35*S5;
    N43XXR5 = -C5*DV35 - S5*WP15;
    N31XXR5 = C4*N41XXR5 + DV25*S4;
    N33XXR5 = -C4*DV25 + N41XXR5*S4;
    N21XXR5 = C3*N31XXR5 + N43XXR5*S3;
    N22XXR5 = -C3*N43XXR5 + N31XXR5*S3;
    N13XXR5 = -C2*N22XXR5 + N21XXR5*S2;
    DG1XXR5 = N13XXR5;
    DG2XXR5 = N33XXR5;
    DG3XXR5 = N33XXR5;
    DG4XXR5 = N43XXR5;
    DG5XXR5 = -DV25;
    No51XY5 = -W15*W35 + WP25;
    No52XY5 = W25*W35 + WP15;
    No53XY5 = W15^2 - W25^2;
    N41XY5 = C5*No51XY5 - No52XY5*S5;
    N43XY5 = -C5*No52XY5 - No51XY5*S5;
    N31XY5 = C4*N41XY5 - No53XY5*S4;
    N33XY5 = C4*No53XY5 + N41XY5*S4;
    N21XY5 = C3*N31XY5 + N43XY5*S3;
    N22XY5 = -C3*N43XY5 + N31XY5*S3;    
    N13XY5 = -C2*N22XY5 + N21XY5*S2;
    DG1XY5 = N13XY5;
    DG2XY5 = N33XY5;
    DG3XY5 = N33XY5;
    DG4XY5 = N43XY5;
    DG5XY5 = No53XY5;
    No51XZ5 = W15*W25 + WP35;
    No52XZ5 = -W15^2 + W35^2;
    No53XZ5 = -W25*W35 + WP15;
    N41XZ5 = C5*No51XZ5 - No52XZ5*S5;
    N43XZ5 = -C5*No52XZ5 - No51XZ5*S5;
    N31XZ5 = C4*N41XZ5 - No53XZ5*S4;
    N33XZ5 = C4*No53XZ5 + N41XZ5*S4;
    N21XZ5 = C3*N31XZ5 + N43XZ5*S3;
    N22XZ5 = -C3*N43XZ5 + N31XZ5*S3;
    N13XZ5 = -C2*N22XZ5 + N21XZ5*S2;
    DG1XZ5 = N13XZ5;
    DG2XZ5 = N33XZ5;
    DG3XZ5 = N33XZ5;
    DG4XZ5 = N43XZ5;
    DG5XZ5 = No53XZ5;
    No51YZ5 = W25^2 - W35^2;
    No52YZ5 = -W15*W25 + WP35;
    No53YZ5 = W15*W35 + WP25;
    N41YZ5 = C5*No51YZ5 - No52YZ5*S5;
    N43YZ5 = -C5*No52YZ5 - No51YZ5*S5;
    N31YZ5 = C4*N41YZ5 - No53YZ5*S4;
    N33YZ5 = C4*No53YZ5 + N41YZ5*S4;
    N21YZ5 = C3*N31YZ5 + N43YZ5*S3;
    N22YZ5 = -C3*N43YZ5 + N31YZ5*S3;
    N13YZ5 = -C2*N22YZ5 + N21YZ5*S2;
    DG1YZ5 = N13YZ5;
    DG2YZ5 = N33YZ5;
    DG3YZ5 = N33YZ5;
    DG4YZ5 = N43YZ5;
    DG5YZ5 = No53YZ5;
    N41ZZR5 = C5*DV55 + DV35*S5;
    N43ZZR5 = C5*DV35 - DV55*S5;
    N31ZZR5 = C4*N41ZZR5 - S4*WP35;
    N33ZZR5 = C4*WP35 + N41ZZR5*S4;
    N21ZZR5 = C3*N31ZZR5 + N43ZZR5*S3;
    N22ZZR5 = -C3*N43ZZR5 + N31ZZR5*S3;
    N13ZZR5 = -C2*N22ZZR5 + N21ZZR5*S2;
    DG1ZZR5 = N13ZZR5;
    DG2ZZR5 = N33ZZR5;
    DG3ZZR5 = N33ZZR5;
    DG4ZZR5 = N43ZZR5;
    DG5ZZR5 = WP35;
    FDI51MX5 = C5*U115 - S5*U215;
    FDI53MX5 = -C5*U215 - S5*U115;
    N41MX5 = S5*VP24;
    N43MX5 = C5*VP24;
    FDI41MX5 = C4*FDI51MX5 - S4*U315;
    FDI43MX5 = C4*U315 + FDI51MX5*S4;
    N31MX5 = C4*N41MX5 - FDI43MX5*RL4 - S4*VP25;
    N32MX5 = -D4*FDI43MX5 - N43MX5;
    N33MX5 = C4*VP25 - D4*FDI53MX5 + FDI41MX5*RL4 + N41MX5*S4;
    FDI32MX5 = -C3*FDI53MX5 + FDI41MX5*S3;
    N21MX5 = C3*N31MX5 - N32MX5*S3;
    N22MX5 = C3*N32MX5 - D3*FDI43MX5 + N31MX5*S3;
    N23MX5 = D3*FDI32MX5 + N33MX5;
    N13MX5 = -C2*N22MX5 - D2*FDI43MX5 + N21MX5*S2;
    DG1MX5 = N13MX5;
    DG2MX5 = N23MX5;
    DG3MX5 = N33MX5;
    DG4MX5 = N43MX5;
    DG5MX5 = VP25;
    FDI51MYR5 = C5*U125 - S5*U225;
    FDI53MYR5 = -C5*U225 - S5*U125;
    FDI41MYR5 = C4*FDI51MYR5 - S4*U325;
    FDI43MYR5 = C4*U325 + FDI51MYR5*S4;
    N31MYR5 = C4*N43MX5 - FDI43MYR5*RL4 + S4*VP15;
    N32MYR5 = -D4*FDI43MYR5 + N41MX5;
    N33MYR5 = -C4*VP15 - D4*FDI53MYR5 + FDI41MYR5*RL4 + N43MX5*S4;
    FDI32MYR5 = -C3*FDI53MYR5 + FDI41MYR5*S3;
    N21MYR5 = C3*N31MYR5 - N32MYR5*S3;
    N22MYR5 = C3*N32MYR5 - D3*FDI43MYR5 + N31MYR5*S3;
    N23MYR5 = D3*FDI32MYR5 + N33MYR5;
    N13MYR5 = -C2*N22MYR5 - D2*FDI43MYR5 + N21MYR5*S2;
    DG1MYR5 = N13MYR5;
    DG2MYR5 = N23MYR5;
    DG3MYR5 = N33MYR5;
    DG4MYR5 = -N41MX5;
    DG5MYR5 = -VP15;
    DG5IA5 = K5^2*QDP5;%*
    % DG5FS5 = sign(QP5);
    % DG5FV5 = QP5;
    N51XXR6 = C6*WP16 - DV36*S6;
    N53XXR6 = C6*DV36 + S6*WP16;
    N41XXR6 = C5*N51XXR6 - DV26*S5;
    N43XXR6 = -C5*DV26 - N51XXR6*S5;
    N31XXR6 = C4*N41XXR6 - N53XXR6*S4;
    N33XXR6 = C4*N53XXR6 + N41XXR6*S4;
    N21XXR6 = C3*N31XXR6 + N43XXR6*S3;
    N22XXR6 = -C3*N43XXR6 + N31XXR6*S3;
    N13XXR6 = -C2*N22XXR6 + N21XXR6*S2;
    DG1XXR6 = N13XXR6;
    DG2XXR6 = N33XXR6;
    DG3XXR6 = N33XXR6;
    DG4XXR6 = N43XXR6;
    DG5XXR6 = N53XXR6;
    DG6XXR6 = -DV26;
    No61XY6 = -W16*W36 + WP26;
    No62XY6 = W26*W36 + WP16;
    No63XY6 = W16^2 - W26^2;
    N51XY6 = C6*No61XY6 - No62XY6*S6;
    N53XY6 = C6*No62XY6 + No61XY6*S6;
    N41XY6 = C5*N51XY6 + No63XY6*S5;
    N43XY6 = C5*No63XY6 - N51XY6*S5;
    N31XY6 = C4*N41XY6 - N53XY6*S4;
    N33XY6 = C4*N53XY6 + N41XY6*S4;
    N21XY6 = C3*N31XY6 + N43XY6*S3;
    N22XY6 = -C3*N43XY6 + N31XY6*S3;
    N13XY6 = -C2*N22XY6 + N21XY6*S2;
    DG1XY6 = N13XY6;
    DG2XY6 = N33XY6;
    DG3XY6 = N33XY6;
    DG4XY6 = N43XY6;
    DG5XY6 = N53XY6;
    DG6XY6 = No63XY6;
    No61XZ6 = W16*W26 + WP36;
    No62XZ6 = -W16^2 + W36^2;
    No63XZ6 = -W26*W36 + WP16;
    N51XZ6 = C6*No61XZ6 - No62XZ6*S6;
    N53XZ6 = C6*No62XZ6 + No61XZ6*S6;
    N41XZ6 = C5*N51XZ6 + No63XZ6*S5;
    N43XZ6 = C5*No63XZ6 - N51XZ6*S5;
    N31XZ6 = C4*N41XZ6 - N53XZ6*S4;
    N33XZ6 = C4*N53XZ6 + N41XZ6*S4;
    N21XZ6 = C3*N31XZ6 + N43XZ6*S3;
    N22XZ6 = -C3*N43XZ6 + N31XZ6*S3;
    N13XZ6 = -C2*N22XZ6 + N21XZ6*S2;
    DG1XZ6 = N13XZ6;
    DG2XZ6 = N33XZ6;
    DG3XZ6 = N33XZ6;
    DG4XZ6 = N43XZ6;
    DG5XZ6 = N53XZ6;
    DG6XZ6 = No63XZ6;
    No61YZ6 = W26^2 - W36^2;
    No62YZ6 = -W16*W26 + WP36;
    No63YZ6 = W16*W36 + WP26;
    N51YZ6 = C6*No61YZ6 - No62YZ6*S6;
    N53YZ6 = C6*No62YZ6 + No61YZ6*S6;
    N41YZ6 = C5*N51YZ6 + No63YZ6*S5;
    N43YZ6 = C5*No63YZ6 - N51YZ6*S5;
    N31YZ6 = C4*N41YZ6 - N53YZ6*S4;
    N33YZ6 = C4*N53YZ6 + N41YZ6*S4;
    N21YZ6 = C3*N31YZ6 + N43YZ6*S3;
    N22YZ6 = -C3*N43YZ6 + N31YZ6*S3;  
    N13YZ6 = -C2*N22YZ6 + N21YZ6*S2;
    DG1YZ6 = N13YZ6;
    DG2YZ6 = N33YZ6;
    DG3YZ6 = N33YZ6;
    DG4YZ6 = N43YZ6;
    DG5YZ6 = N53YZ6;
    DG6YZ6 = No63YZ6;
    N51ZZ6 = C6*DV56 + DV36*S6;
    N53ZZ6 = -C6*DV36 + DV56*S6;
    N41ZZ6 = C5*N51ZZ6 + S5*WP36;
    N43ZZ6 = C5*WP36 - N51ZZ6*S5;
    N31ZZ6 = C4*N41ZZ6 - N53ZZ6*S4;
    N33ZZ6 = C4*N53ZZ6 + N41ZZ6*S4;
    N21ZZ6 = C3*N31ZZ6 + N43ZZ6*S3;
    N22ZZ6 = -C3*N43ZZ6 + N31ZZ6*S3;
    N13ZZ6 = -C2*N22ZZ6 + N21ZZ6*S2;
    DG1ZZ6 = N13ZZ6;
    DG2ZZ6 = N33ZZ6;
    DG3ZZ6 = N33ZZ6;
    DG4ZZ6 = N43ZZ6;
    DG5ZZ6 = N53ZZ6;
    DG6ZZ6 = WP36;
    FDI61MX6 = C6*U116 - S6*U216;
    FDI63MX6 = C6*U216 + S6*U116;
    N51MX6 = -FDI63MX6*RL6 - S6*VSP26;
    N53MX6 = C6*VSP26 + FDI61MX6*RL6;
    FDI51MX6 = C5*FDI61MX6 + S5*U316;
    FDI53MX6 = C5*U316 - FDI61MX6*S5;
    N41MX6 = C5*N51MX6 + S5*VP26;
    N43MX6 = C5*VP26 - N51MX6*S5;
    FDI41MX6 = C4*FDI51MX6 - FDI63MX6*S4;
    FDI43MX6 = C4*FDI63MX6 + FDI51MX6*S4;
    N31MX6 = C4*N41MX6 - FDI43MX6*RL4 - N53MX6*S4;
    N32MX6 = -D4*FDI43MX6 - N43MX6;
    N33MX6 = C4*N53MX6 - D4*FDI53MX6 + FDI41MX6*RL4 + N41MX6*S4;
    FDI32MX6 = -C3*FDI53MX6 + FDI41MX6*S3;
    N21MX6 = C3*N31MX6 - N32MX6*S3;
    N22MX6 = C3*N32MX6 - D3*FDI43MX6 + N31MX6*S3;
    N23MX6 = D3*FDI32MX6 + N33MX6;
    N13MX6 = -C2*N22MX6 - D2*FDI43MX6 + N21MX6*S2;
    DG1MX6 = N13MX6;
    DG2MX6 = N23MX6;
    DG3MX6 = N33MX6;
    DG4MX6 = N43MX6;
    DG5MX6 = N53MX6;
    DG6MX6 = VP26;
    FDI61MY6 = C6*U126 - S6*U226;
    FDI63MY6 = C6*U226 + S6*U126;
    N51MY6 = -C6*VSP26 - FDI63MY6*RL6;
    N53MY6 = FDI61MY6*RL6 - S6*VSP26;
    FDI51MY6 = C5*FDI61MY6 + S5*U326;
    FDI53MY6 = C5*U326 - FDI61MY6*S5;
    N41MY6 = C5*N51MY6 - S5*VP16;
    N43MY6 = -C5*VP16 - N51MY6*S5;
    FDI41MY6 = C4*FDI51MY6 - FDI63MY6*S4;
    FDI43MY6 = C4*FDI63MY6 + FDI51MY6*S4;
    N31MY6 = C4*N41MY6 - FDI43MY6*RL4 - N53MY6*S4;
    N32MY6 = -D4*FDI43MY6 - N43MY6;
    N33MY6 = C4*N53MY6 - D4*FDI53MY6 + FDI41MY6*RL4 + N41MY6*S4;
    FDI32MY6 = -C3*FDI53MY6 + FDI41MY6*S3;
    N21MY6 = C3*N31MY6 - N32MY6*S3;
    N22MY6 = C3*N32MY6 - D3*FDI43MY6 + N31MY6*S3;
    N23MY6 = D3*FDI32MY6 + N33MY6;
    N13MY6 = -C2*N22MY6 - D2*FDI43MY6 + N21MY6*S2;
    DG1MY6 = N13MY6;
    DG2MY6 = N23MY6;
    DG3MY6 = N33MY6;
    DG4MY6 = N43MY6;
    DG5MY6 = N53MY6;
    DG6MY6 = -VP16;

    % Ratio56=1.116279;
    % Ratio56=0;
    % DG6IA6 = QDP6 + Ratio56*QDP5;  %耦合
    % DG6FS6 = sign(QP6);%耦合
    % DG6FV6 = QP6;

    % DG6IA6 = K6^2*(QDP6+QDP5); %*
    DG6IA6 = K6^2*(QDP6); %*
    % ee=0.03;
 
     if abs(QP6)<0.7 
        DG6FS61 = sign(QP6 );
        DG6FV61 = QP6 ;

        DG6FS62 = 0;
        DG6FV62 = 0 ;
        % DG6FS63 = 0;
        % DG6FV63 = 0 ;        
    elseif abs(QP6)>=0.7
        DG6FS61 = 0;
        DG6FV61 = 0 ;

        DG6FS62 = sign(QP6 );
        DG6FV62 = QP6 ;
        % DG6FS63 = 0;
        % DG6FV63 = 0 ; 
    end
    % % 摩擦正负相同 56参数
    if abs(QP5)<0.7 
        DG5FS51 = sign(QP5 );
        DG5FV51 = QP5 ;

        DG5FS52 = 0;
        DG5FV52 = 0 ;
    elseif abs(QP5)>=0.7
        DG5FS51 = 0;
        DG5FV51 = 0 ;

        DG5FS52 = sign(QP5 );
        DG5FV52 = QP5 ;
    end
   % % 摩擦正负相同 56参数
   if abs(QP4)<0.7 
        DG4FS41 = sign(QP4 );
        DG4FV41 = QP4 ;

        DG4FS42 = 0;
        DG4FV42 = 0 ;
        % DG4FS43 = 0;
        % DG4FV43 = 0 ;
    elseif abs(QP4)>=0.7
        DG4FS41 = 0;
        DG4FV41 = 0 ;

        DG4FS42 = sign(QP4 );
        DG4FV42 = QP4 ;
        % DG4FS43 = 0;
        % DG4FV43 = 0 ;        
    end
    % % 摩擦正负相同 56参数
    % if abs(QP3)<0.7 
    %     DG3FS31 = sign(QP3 );
    %     DG3FV31 = QP3 ;
    % 
    %     DG3FS32 = 0;
    %     DG3FV32 = 0 ;
    % elseif abs(QP3)>=0.7
    %     DG3FS31 = 0;
    %     DG3FV31 = 0 ;
    % 
    %     DG3FS32 = sign(QP3 );
    %     DG3FV32 = QP3 ;
    % end
    % % 摩擦正负相同 56参数
    % if abs(QP2)<0.7 
    %     DG2FS21 = sign(QP2 );
    %     DG2FV21 = QP2 ;
    % 
    %     DG2FS22 = 0;
    %     DG2FV22 = 0 ;
    % elseif abs(QP2)>=0.7
    %     DG2FS21 = 0;
    %     DG2FV21 = 0 ;
    % 
    %     DG2FS22 = sign(QP2 );
    %     DG2FV22 = QP2 ;
    % end    
    % % 摩擦正负相同 56参数
    if abs(QP1)<0.7 
        DG1FS11 = sign(QP1 );
        DG1FV11 = QP1 ;

        DG1FS12 = 0;
        DG1FV12 = 0 ;
    elseif abs(QP1)>=0.7
        DG1FS11 = 0;
        DG1FV11 = 0 ;

        DG1FS12 = sign(QP1 );
        DG1FV12 = QP1 ;
    end  

    w_min11=DG1ZZR1;            %ZZR1
    w_min12=DG1XXR2;            %XXR2
    w_min13=DG1XY2;             %XY2
    w_min14=DG1XZR2;            %XZR2
    w_min15=DG1YZ2;             %YZ2
    w_min16=DG1ZZR2;            %ZZR2
    w_min17=DG1MXR2;            %MXR2
    w_min18=DG1MY2;             %MY2
    w_min19=DG1XXR3;            %XXR3
    w_min110=DG1XYR3;           %XYR3
    w_min111=DG1XZ3;            %XZ3
    w_min112=DG1YZ3;            %YZ3
    w_min113=DG1ZZR3;           %ZZR3
    w_min114=DG1MXR3;           %MXR3
    w_min115=DG1MYR3;           %MYR3
    w_min116=0;     % DG1IA3     %IA3
    w_min117=DG1XXR4;          %XXR4
    w_min118=DG1XY4;            %XY4
    w_min119=DG1XZ4;            %XZ4
    w_min120=DG1YZ4;            %YZ4
    w_min121=DG1ZZR4;           %ZZR4
    w_min122=DG1MX4;            %MX4
    w_min123=DG1MYR4;           %MYR4
    w_min124=0;        %DG1IA4     %IA4
    w_min125=DG1XXR5;           %XXR5
    w_min126=DG1XY5;            %XY5
    w_min127=DG1XZ5;            %XZ5
    w_min128=DG1YZ5;            %YZ5
    w_min129=DG1ZZR5;           %ZZR5
    w_min130=DG1MX5;            %MX5
    w_min131=DG1MYR5;           %MYR5
    w_min132=0;         %DG1IA5     %IA5
    w_min133=DG1XXR6;           %XXR6
    w_min134=DG1XY6;            %XY6
    w_min135=DG1XZ6;            %XZ6
    w_min136=DG1YZ6;            %YZ6
    w_min137=DG1ZZ6;            %ZZ6
    w_min138=DG1MX6;            %MX6
    w_min139=DG1MY6;            %MY6
    w_min140=0;                 %DG1IA6     %IA6
    w_min141=0;                 %FS56
    w_min142=0;          %FV56    
    w_min143=DG1FS11;            %FS1
    w_min144=DG1FV11;            %FV1
    w_min145=DG1FS12;   %DG1FS2         %FS2
    w_min146=DG1FV12;         %DG1FV2         %FV2
    w_min147=0;         %DG1FS3         %FS3
    w_min148=0;         %DG1FV3         %FV3
    w_min149=0;         %DG1FS4         %FS4
    w_min150=0;         %DG1FV4         %FV4
    w_min151=0;         %DG1FS5         %FS5
    w_min152=0;         %DG1FV5         %FV5
    w_min153=0;         %DG1FS6         %FS6
    w_min154=0;           %DG1FV6         %FV6                     
    w_min155=0;                
    w_min156=0;
    w_min157=0;                
    w_min158=0;
    w_min159=0;                
    w_min160=0;    
    w_min161=0;                
    w_min162=0;    
    % w_min163=0;                
    % w_min164=0;  
    % w_min165=0;                
    % w_min166=0; 
    % w_min153=1;                 %F01   grotjahn friction
    % w_min154=0;                 %F02
    % w_min155=0;                 %F03
    % w_min156=0;                 %F04
    % w_min157=0;                 %F05
    % w_min158=0;                 %F06
    % 
    % w_min159=QP1^(1/3);         %F31
    % w_min160=0;                 %F32
    % w_min161=0;                 %F33
    % w_min162=0;                 %F34
    % w_min163=0;                 %F35
    % w_min164=0;                 %F36

    w_min21=0; %DG2ZZR1
    w_min22=DG2XXR2;
    w_min23=DG2XY2;
    w_min24=DG2XZR2;
    w_min25=DG2YZ2;
    w_min26=DG2ZZR2;
    w_min27=DG2MXR2;
    w_min28=DG2MY2;
    w_min29=DG2XXR3;
    w_min210=DG2XYR3;
    w_min211=DG2XZ3;
    w_min212=DG2YZ3;
    w_min213=DG2ZZR3;
    w_min214=DG2MXR3;
    w_min215=DG2MYR3;
    w_min216=0;   %DG2IA3
    w_min217=DG2XXR4;
    w_min218=DG2XY4;
    w_min219=DG2XZ4;
    w_min220=DG2YZ4;
    w_min221=DG2ZZR4;
    w_min222=DG2MX4;
    w_min223=DG2MYR4;
    w_min224=0;    %DG2IA4
    w_min225=DG2XXR5;
    w_min226=DG2XY5;
    w_min227=DG2XZ5;
    w_min228=DG2YZ5;
    w_min229=DG2ZZR5;
    w_min230=DG2MX5;
    w_min231=DG2MYR5;
    w_min232=0;     %DG2IA5
    w_min233=DG2XXR6;
    w_min234=DG2XY6;
    w_min235=DG2XZ6;
    w_min236=DG2YZ6;
    w_min237=DG2ZZ6;
    w_min238=DG2MX6;
    w_min239=DG2MY6;
    w_min240=0;     %DG2IA6
    w_min241=0;
    w_min242=0;

    w_min243=0; %DG2FS2
    w_min244=0; %DG2FV2
    w_min245=0; %DG2FS3
    w_min246=0; %DG2FV3
    w_min247=DG2FS2; %DG2FS4
    w_min248=DG2FV2; %DG2FV4
    w_min249=0; %DG2FS5
    w_min250=0; %DG2FV5
    w_min251=0; %DG2FS6
    w_min252=0; %DG2FV6
    w_min253=0;              
    w_min254=0;  
    w_min255=0;
    w_min256=0;                
    w_min257=0;   
    w_min258=0;                
    w_min259=0;  
    w_min260=0;                
    w_min261=0;
    w_min262=0;                
    % w_min263=0;                
    % w_min264=0;  
    % w_min265=0;                
    % w_min266=0;      
    % w_min253=0;                 %F01   grotjahn friction
    % w_min254=1;                 %F02
    % w_min255=0;                 %F03
    % w_min256=0;                 %F04
    % w_min257=0;                 %F05
    % w_min258=0;                 %F06
    % 
    % w_min259=0;                 %F31
    % w_min260=QP2^(1/3);         %F32
    % w_min261=0;                 %F33
    % w_min262=0;                 %F34
    % w_min263=0;                 %F35
    % w_min264=0;                 %F36


    w_min31=0; %DG3ZZR1
    w_min32=0; %DG3XXR2
    w_min33=0; %DG3XY2
    w_min34=0; %DG3XZR2
    w_min35=0; %DG3YZ2
    w_min36=0; %DG3ZZR2
    w_min37=0; %DG3MXR2
    w_min38=0; %DG3MXR2
    w_min39=DG3XXR3;
    w_min310=DG3XYR3;
    w_min311=DG3XZ3;
    w_min312=DG3YZ3;
    w_min313=DG3ZZR3;
    w_min314=DG3MXR3;
    w_min315=DG3MYR3;
    w_min316=DG3IA3;  

    w_min317=DG3XXR4;
    w_min318=DG3XY4;
    w_min319=DG3XZ4;
    w_min320=DG3YZ4;
    w_min321=DG3ZZR4;
    w_min322=DG3MX4;
    w_min323=DG3MYR4;
    w_min324=0;    %DG3IA4
    w_min325=DG3XXR5;
    w_min326=DG3XY5;
    w_min327=DG3XZ5;
    w_min328=DG3YZ5;
    w_min329=DG3ZZR5;
    w_min330=DG3MX5;
    w_min331=DG3MYR5;
    w_min332=0;     %DG3IA5
    w_min333=DG3XXR6;
    w_min334=DG3XY6;
    w_min335=DG3XZ6;
    w_min336=DG3YZ6;
    w_min337=DG3ZZ6;
    w_min338=DG3MX6;
    w_min339=DG3MY6;
    w_min340=0;     %DG3IA6
    w_min341=0;     %DG3FS56    
    w_min342=0;     %DG3Fv56
    w_min343=0; %DG3FS1
    w_min344=0; %DG3FV1
    w_min345=0; %DG3FS3
    w_min346=0; %DG3FV3
    w_min347=0; %DG3FS4
    w_min348=0; %DG3FV4
    w_min349=DG3FS3; %DG3FS5
    w_min350=DG3FV3; %DG3FV5
    w_min351=0; %DG3FS6
    w_min352=0; %DG3FV6
    w_min353=0;                
    w_min354=0;
    w_min355=0;
    w_min356=0;                
    w_min357=0;
    w_min358=0;                
    w_min359=0;   
    w_min360=0;                
    w_min361=0;
    w_min362=0;                
    % w_min363=0; 
    % w_min364=0;
    % w_min365=0;
    % w_min366=0;
    % w_min353=0;                 %F01   grotjahn friction
    % w_min354=0;                 %F02
    % w_min355=1;                 %F03
    % w_min356=0;                 %F04
    % w_min357=0;                 %F05
    % w_min358=0;                 %F06
    % 
    % w_min359=0;                 %F31
    % w_min360=0;                 %F32
    % w_min361=QP3^(1/3);         %F33
    % w_min362=0;                 %F34
    % w_min363=0;                 %F35
    % w_min364=0;                 %F36

    w_min41=0; %DG4ZZR1
    w_min42=0; %DG4XXR2
    w_min43=0; %DG4XY2
    w_min44=0; %DG4XZR2
    w_min45=0; %DG4YZ2
    w_min46=0; %DG4ZZR2
    w_min47=0; %DG4MXR2
    w_min48=0; %DG4MXR2
    w_min49=0; %DG4XXR3
    w_min410=0; %DG4XYR3
    w_min411=0; %DG4XZ3
    w_min412=0; %DG4YZ3
    w_min413=0; %DG4ZZR3
    w_min414=0; %DG4MXR3
    w_min415=0; %DG4MYR3
    w_min416=0;
    w_min417=DG4XXR4;%DG4IA3  
    w_min418=DG4XY4;
    w_min419=DG4XZ4;
    w_min420=DG4YZ4;
    w_min421=DG4ZZR4;
    w_min422=DG4MX4;
    w_min423=DG4MYR4;
    w_min424=DG4IA4;    %DG4IA4
    w_min425=DG4XXR5;
    w_min426=DG4XY5;
    w_min427=DG4XZ5;
    w_min428=DG4YZ5;
    w_min429=DG4ZZR5;
    w_min430=DG4MX5;
    w_min431=DG4MYR5;
    w_min432=0;     %DG4IA5
    w_min433=DG4XXR6;
    w_min434=DG4XY6;
    w_min435=DG4XZ6;
    w_min436=DG4YZ6;
    w_min437=DG4ZZ6;
    w_min438=DG4MX6;
    w_min439=DG4MY6;
    w_min440=0;     %DG4IA6    %DG4FS56    
    w_min441=0;     %DG4Fv56
    w_min442=0; %DG4FS2
    w_min443=0; %DG4FV2
    w_min444=0; %DG4FS3
    w_min445=0; %DG4FV3
    w_min446=0; %COUPLE56
    w_min447=0; %COUPLE56
    w_min448=0; %DG4FS5
    w_min449=0; %DG4FV5
    w_min450=0;
    w_min451=DG4FS41;
    w_min452=DG4FV41;
    w_min453=DG4FS42;
    w_min454=DG4FV42;
    w_min455=0;                 
    w_min456=0;      
    w_min457=0;   
    w_min458=0;      
    w_min459=0;            
    w_min460=0;      
    w_min461=0;   
    w_min462=0; 
    % w_min463=0; 
    % w_min464=0; 
    % w_min465=0; 
    % w_min466=0; 

    % w_min453=0;                 %F01   grotjahn friction
    % w_min454=0;                 %F02
    % w_min455=0;                 %F03
    % w_min456=1;                 %F04
    % w_min457=0;                 %F05
    % w_min458=0;                 %F06
    % 
    % w_min459=0;                 %F31
    % w_min460=0;                 %F32
    % w_min461=0;                 %F33
    % w_min462=QP4^(1/3);         %F34
    % w_min463=0;                 %F35
    % w_min464=0;                 %F36

    w_min51=0; %DG5ZZR1
    w_min52=0; %DG5XXR2
    w_min53=0; %DG5XY2
    w_min54=0; %DG5XZR2
    w_min55=0; %DG5YZ2
    w_min56=0; %DG5ZZR2
    w_min57=0; %DG5MXR2
    w_min58=0; %DG5MXR2
    w_min59=0; %DG5XXR3
    w_min510=0; %DG5XYR3
    w_min511=0; %DG5XZ3
    w_min512=0; %DG5YZ3
    w_min513=0; %DG5ZZR3
    w_min514=0; %DG5MXR3
    w_min515=0; %DG5MYR3
    w_min516=0; %DG5XXR4
    w_min517=0; %DG5XY4
    w_min518=0; %DG5XZ4
    w_min519=0; %DG5YZ4
    w_min520=0; %DG5ZZR4
    w_min521=0; %DG5MX4
    w_min522=0; %DG5MYR4
    w_min523=0;    %DG5IA4
    w_min524=0;
    w_min525=DG5XXR5;
    w_min526=DG5XY5;
    w_min527=DG5XZ5;
    w_min528=DG5YZ5;
    w_min529=DG5ZZR5;
    w_min530=DG5MX5;
    w_min531=DG5MYR5;
    w_min532=DG5IA5;    
    w_min533=DG5XXR6;
    w_min534=DG5XY6;
    w_min535=DG5XZ6;
    w_min536=DG5YZ6;
    w_min537=DG5ZZ6;
    w_min538=DG5MX6;
    w_min539=DG5MY6;
    w_min540=0;     %DG5IA6  **
    % w_min540=(K6)^2*(QDP6+QDP5); 
    % w_min540=(K6)^2*(QDP6); 
    w_min541=sign(QP6);       %DG5FS1    
    w_min542=QP6;             %DG5Fv1 

                 %DG5FS2
    w_min543=0; %DG5FV2
    w_min544=0; %DG5FS3
    w_min545=0; %DG5FV3
    w_min546=0; %DG5FS4

    w_min547=0; %DG5FV4
    w_min548=0; %DG5FS5
    w_min549=0; %DG5FV5
    w_min550=0; %DG5FS6

    w_min551=0; %DG5FV6  
    w_min552=0;              
    w_min553=0;
    w_min554=0;
    w_min555=DG5FS51;
    w_min556=DG5FV51;          
    w_min557=DG5FS52;
    w_min558=DG5FV52;             
    w_min559=0;
    w_min560=0;    
    w_min561=0;
    w_min562=0;   

    % w_min563=0;   
    % w_min564=0; 
    % w_min565=0; 
    % w_min566=0; 


    % w_min553=0;                 %F01   grotjahn friction
    % w_min554=0;                 %F02
    % w_min555=0;                 %F03
    % w_min556=0;                 %F04
    % w_min557=1;                 %F05
    % w_min558=0;                 %F06
    % 
    % w_min559=0;                 %F31
    % w_min560=0;                 %F32
    % w_min561=0;                 %F33
    % w_min562=0;                 %F34
    % w_min563=QP5^(1/3);         %F35
    % w_min564=0;                 %F36

    w_min61=0; %DG6ZZR1
    w_min62=0; %DG6XXR2
    w_min63=0; %DG6XY2
    w_min64=0; %DG6XZR2
    w_min65=0; %DG6YZ2
    w_min66=0; %DG6ZZR2
    w_min67=0; %DG6MXR2
    w_min68=0; %DG6MXR2
    w_min69=0; %DG6XXR3
    w_min610=0; %DG6XYR3
    w_min611=0; %DG6XZ3
    w_min612=0; %DG6YZ3
    w_min613=0; %DG6ZZR3
    w_min614=0; %DG6MXR3
    w_min615=0; %DG6MYR3
    w_min616=0; %DG6XXR4 %DG6IA3  
    w_min617=0; %DG6XY4
    w_min618=0; %DG6XZ4
    w_min619=0; %DG6YZ4
    w_min620=0; %DG6ZZR4
    w_min621=0; %DG6MX4
    w_min622=0; %DG6MYR4
    w_min623=0;    %DG6IA4
    w_min624=0; %DG6XXR5
    w_min625=0; %DG6XY5
    w_min626=0; %DG6XZ5
    w_min627=0; %DG6YZ5
    w_min628=0; %DG6ZZR5
    w_min629=0; %DG6MX5
    w_min630=0; %DG6MYR5
    w_min631=0;     %DG6IA5
    w_min632=0;
    w_min633=DG6XXR6;
    w_min634=DG6XY6;
    w_min635=DG6XZ6;
    w_min636=DG6YZ6;
    w_min637=DG6ZZ6;
    w_min638=DG6MX6;
    w_min639=DG6MY6; %DG6IA6
    w_min640=DG6IA6;       %DG6FS56    
    w_min641=0;     %DG6Fv56
    w_min642=0; %DG6FS2
    w_min643=0; %DG6FV2
    w_min644=0; %DG6FS3
    w_min645=0; %DG6FV3
    w_min646=0; %DG6FS4
    w_min647=0; %DG6FV4
    w_min648=0; %DG6FS5
    w_min649=0; %DG6FV5
    w_min650=0;
    w_min651=0;
    w_min652=0;
    w_min653=0;
    w_min654=0;
    w_min655=0;
    w_min656=0;
    w_min657=0;
    w_min658=0;
    w_min659=DG6FS61; 
    w_min660=DG6FV61;
    w_min661=DG6FS62; 
    w_min662=DG6FV62;

    % w_min663=0;
    % w_min664=0;  
    % w_min665=DG6FS63;
    % w_min666=DG6FV63;

    % % w_min653= -sign(QP6 + 1.116279*QP5);       %FS56
    % % w_min654= -QP6 - 1.116279*QP5;   
    % % w_min653=sign(dq1(6));         %FS56
    % % w_min654=dq1(6);
    % w_min653=DG6FV62;         %FS56
    % w_min654=DG6FS63;
    % w_min655=DG6FV63;               
    % w_min656=DG6FS64;  
    % w_min657=DG6FS65;            
    % w_min658=DG6FV65;  
    % w_min659=DG6FS66;             
    % w_min660=DG6FV66;  
    % w_min661=0;             
    % w_min662=0;  



    % w_min653=sign(QP6 - 1.116279*QP5 );                 %F01   grotjahn friction
    % w_min654=(QP6 - 1.116279*QP5 );                 %F02
    % w_min655=0;                 %F03
    % w_min656=0;                 %F04
    % w_min657=0;                 %F05
    % w_min658=1;                 %F06
    % 
    % w_min659=0;                 %F31
    % w_min660=0;                 %F32
    % w_min661=0;                 %F33
    % w_min662=0;                 %F34
    % w_min663=0;                 %F35
    % w_min664=QP6^(1/3);         %F36


    w_min(1, :) = [w_min11, w_min12, w_min13, w_min14, w_min15, w_min16, w_min17, w_min18, w_min19, w_min110, ...
                   w_min111, w_min112, w_min113, w_min114, w_min115, w_min116, w_min117, w_min118, w_min119, w_min120, ...
			       w_min121, w_min122, w_min123, w_min124, w_min125, w_min126, w_min127, w_min128, w_min129, w_min130, ...
			       w_min131, w_min132, w_min133, w_min134, w_min135, w_min136, w_min137, w_min138, w_min139, w_min140, ...
			       w_min141, w_min142,w_min143, w_min144, w_min145, w_min146, w_min147, w_min148, w_min149, w_min150, ...
			       w_min151, w_min152, w_min153, w_min154, w_min155, w_min156,w_min157, w_min158, w_min159, w_min160, ...
                   w_min161, w_min162];%,w_min163, w_min164, w_min165, w_min166];                
                   
    w_min(2, :) = [w_min21, w_min22, w_min23, w_min24, w_min25, w_min26, w_min27, w_min28, w_min29, w_min210, ...
                   w_min211, w_min212, w_min213, w_min214, w_min215, w_min216, w_min217, w_min218, w_min219, w_min220, ...
			       w_min221, w_min222, w_min223, w_min224, w_min225, w_min226, w_min227, w_min228, w_min229, w_min230, ...
			       w_min231, w_min232, w_min233, w_min234, w_min235, w_min236, w_min237, w_min238, w_min239, w_min240, ...
			       w_min241, w_min242, w_min243, w_min244, w_min245, w_min246, w_min247, w_min248, w_min249, w_min250, ...
			       w_min251, w_min252, w_min253, w_min254, w_min255, w_min256,w_min257, w_min258, w_min259, w_min260, ...
                   w_min261, w_min262];%,w_min263, w_min264, w_min265, w_min266];  

    w_min(3, :) = [w_min31, w_min32, w_min33, w_min34, w_min35, w_min36, w_min37, w_min38, w_min39, w_min310, ...
                   w_min311, w_min312, w_min313, w_min314, w_min315, w_min316, w_min317, w_min318, w_min319, w_min320, ...
			       w_min321, w_min322, w_min323, w_min324, w_min325, w_min326, w_min327, w_min328, w_min329, w_min330, ...
			       w_min331, w_min332, w_min333, w_min334, w_min335, w_min336, w_min337, w_min338, w_min339, w_min340, ...
			       w_min341, w_min342, w_min343, w_min344, w_min345, w_min346, w_min347, w_min348, w_min349, w_min350, ...
			       w_min351, w_min352, w_min353, w_min354, w_min355, w_min356,w_min357, w_min358, w_min359, w_min360, ...
                   w_min361, w_min362];%,w_min363, w_min364, w_min365, w_min366];  

    w_min(4, :) = [w_min41, w_min42, w_min43, w_min44, w_min45, w_min46, w_min47, w_min48, w_min49, w_min410, ...
                   w_min411, w_min412, w_min413, w_min414, w_min415, w_min416, w_min417, w_min418, w_min419, w_min420, ...
			       w_min421, w_min422, w_min423, w_min424, w_min425, w_min426, w_min427, w_min428, w_min429, w_min430, ...
			       w_min431, w_min432, w_min433, w_min434, w_min435, w_min436, w_min437, w_min438, w_min439, w_min440, ...
			       w_min441, w_min442,w_min443, w_min444, w_min445, w_min446, w_min447, w_min448, w_min449, w_min450, ...
			       w_min451, w_min452, w_min453, w_min454, w_min455, w_min456,w_min457, w_min458, w_min459, w_min460, ...
                   w_min461, w_min462];%,w_min463, w_min464, w_min465, w_min466];  

    w_min(5, :) = [w_min51, w_min52, w_min53, w_min54, w_min55, w_min56, w_min57, w_min58, w_min59, w_min510, ...
                   w_min511, w_min512, w_min513, w_min514, w_min515, w_min516, w_min517, w_min518, w_min519, w_min520, ...
			       w_min521, w_min522, w_min523, w_min524, w_min525, w_min526, w_min527, w_min528, w_min529, w_min530, ...
			       w_min531, w_min532, w_min533, w_min534, w_min535, w_min536, w_min537, w_min538, w_min539, w_min540, ...
			       w_min541, w_min542,w_min543, w_min544, w_min545, w_min546, w_min547, w_min548, w_min549, w_min550, ...
			       w_min551, w_min552,w_min553, w_min554, w_min555, w_min556,w_min557, w_min558, w_min559, w_min560, ...
                   w_min561, w_min562];%,w_min563, w_min564, w_min565, w_min566];  

    w_min(6, :) = [w_min61, w_min62, w_min63, w_min64, w_min65, w_min66, w_min67, w_min68, w_min69, w_min610, ...
                   w_min611, w_min612, w_min613, w_min614, w_min615, w_min616, w_min617, w_min618, w_min619, w_min620, ...
			       w_min621, w_min622, w_min623, w_min624, w_min625, w_min626, w_min627, w_min628, w_min629, w_min630, ...
			       w_min631, w_min632, w_min633, w_min634, w_min635, w_min636, w_min637, w_min638, w_min639, w_min640, ...
			       w_min641, w_min642, w_min643, w_min644, w_min645, w_min646, w_min647, w_min648, w_min649, w_min650, ...
			       w_min651, w_min652,w_min653, w_min654, w_min655, w_min656,w_min657, w_min658, w_min659, w_min660, ...
                   w_min661, w_min662];%,w_min663, w_min664, w_min665, w_min666];  
end

