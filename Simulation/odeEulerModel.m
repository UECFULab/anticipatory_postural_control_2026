%% Calculate next state of internal model (with Euler Method)
function xNext = odeEulerModel(x, u)
global ctrlParam 

xNext    = x +diffEqMPCModel_ms(x, u)*ctrlParam.mpcSamp; % Euler-Maruyama Method
end

%% Internal model for MPC
function dz =diffEqMPCModel_ms(z,mv)
global bodyParam ctrlParam
global currentTime 

% Body parameters
m1  = bodyParam.m1;
J1  = bodyParam.J1;
L1  = bodyParam.L1;
h1  = bodyParam.h1;
m2  = bodyParam.m2;
J2  = bodyParam.J2;
h2  = bodyParam.h2;
g   = bodyParam.g;

time   = z(5); % time

%-- Floor motion --%
withFloorMotion = 1;
if(withFloorMotion)
    soundStart= 13;  % (s)
    tiltStart = 15;  % (s)
    tiltEnd   = 16;  % (s)
    if(currentTime>soundStart) % Start cue: 13 s
        tiltAngle = bodyParam.tiltMax/180*pi; % Maximum tilt angle (default: 4 degree)
        if(time==tiltStart)
            phi      = 0;
            dPhi     = tiltAngle/(tiltEnd-tiltStart); % angular velocity
            lastDPhi = 0;
        elseif(time>tiltStart &&time<tiltEnd)
            phi  =tiltAngle/(tiltEnd-tiltStart)*(time-tiltStart);
            dPhi =tiltAngle/(tiltEnd-tiltStart);
            lastDPhi = dPhi;
        elseif(time==tiltEnd)
            phi = tiltAngle;
            dPhi = 0;
            lastDPhi = tiltAngle/(tiltEnd-tiltStart); % angular velocity
        elseif(time>tiltEnd)
            phi    = tiltAngle;
            dPhi   = 0;
            lastDPhi = dPhi;
        else
            phi  = 0;
            dPhi = 0;
            lastDPhi = dPhi;
        end
    else
        phi  = 0;
        dPhi = 0;
        lastDPhi = 0;
    end
else
    phi   = 0;
    dPhi  = 0;
    lastDPhi = 0;
end
z(6) = phi;

% State variables
dqDList = [z(1);z(3)]; % dq1 and dq2
qDList  = [z(2);z(4)]; % q1 and q2
qD1 = qDList(1);
qD2 = qDList(2);
dqD1 = dqDList(1);
dqD2 = dqDList(2);

%========== Setting muscles and torques ==========%
mv1 = mv(1);
mv2 = mv(2);
mv3 = mv(3);
%
% Active muscle forces
[mTA_fm, mTA_lm, mTA_dLm]       = calcMuscleForce( mv1, 'TA',   qD1, dqD1, true); % TA
[mGC_fm, mGC_lm, mGC_dLm]       = calcMuscleForce( mv2, 'GC',   qD1, dqD1, false); % GC
[mIP_fm, mIP_lm, mIP_dLm]       = calcMuscleForce( mv3, 'IP',   qD2, dqD2, true); 
[mGMax_fm, mGMax_lm, mGMax_dLm] = calcMuscleForce(-mv3, 'GMax', qD2, dqD2, false); 
%
% Passive muscle forces due to muscle stiffness control
fmStiffTA   = calcFmStiff(mTA_lm,   mTA_dLm,   'TA', true);
fmStiffGC   = calcFmStiff(mGC_lm,   mGC_dLm,   'GC', true);
fmStiffIP   = calcFmStiff(mIP_lm,   mIP_dLm,   'IP', false);
fmStiffGMax = calcFmStiff(mGMax_lm, mGMax_dLm, 'GMax', false);
%
% Calculation of torque by muscles
%- Active torque -%
uA = -bodyParam.TA_ma*mTA_fm +bodyParam.GC_ma*mGC_fm;
uH = -bodyParam.IP_ma*mIP_fm +bodyParam.GMax_ma*mGMax_fm;
%- Passive torque -%
uAStiff = -bodyParam.TA_ma*fmStiffTA +bodyParam.GC_ma*fmStiffGC;
uHStiff = -bodyParam.IP_ma*fmStiffIP +bodyParam.GMax_ma*fmStiffGMax;
%==================================================%

%==========     Equations of motion      ==========%
%- M term (Inertia)-%
M11 = J1 +J2 +m1*(h1^2) +m2*((L1^2)+(h2^2))+2*m2*L1*h2*cos(qD2);
M12 = J2 +m2*h2^2+m2*L1*h2*cos(qD2);
M21 = J2 +m2*h2^2+m2*L1*h2*cos(qD2);
M22 = J2+m2*h2^2;
mMat = [M11 M12; M21 M22];
%
%- C term (Centrifugal force and Coriolis force) -%
C1_11 = m2*L1*h2*sin(qD2)*(-2*dqD1*dqD2 -(dqD2^2) -2*dqD2*dPhi);
C1_21 = m2*L1*h2*sin(qD2)*(dqD1^2 + 2*dqD1*dPhi +(dPhi^2));
c1Mat = [C1_11; C1_21];
%
%- G term (Gravity) -%
G11 = -((m1*h1+m2*L1)*sin(qD1+phi)+m2*h2*sin(qD1+qD2+phi));
G21 = -m2*h2*sin(qD1+qD2+phi);
gMat= [G11; G21];
%
% Active Torque term -%
ctrlTorque = [uA-uH; uH]; % no noise
% Passive Torque term -%
tMat = [uAStiff-uHStiff; uHStiff];
%
%-- Equation of motion --%
mEqDDQ = pinv(mMat)*(-gMat*g -c1Mat +ctrlTorque +tMat -[1;0]*(dPhi-lastDPhi)/ctrlParam.mpcSamp);
mEqDQ = dqDList;
%
dz = [mEqDDQ(1); mEqDQ(1); mEqDDQ(2); mEqDQ(2); 1; dPhi]; 
end

%% Calculation of muscle force
function [fm, lm, dLm] = calcMuscleForce(act, muscleName, qD, dqD, isFront)
global bodyParam

fmBar = bodyParam.([muscleName '_fmBar']);
lmBar = bodyParam.([muscleName '_lmBar']);
dLmBar= bodyParam.([muscleName '_dLmBar']);
ma= bodyParam.([muscleName '_ma']);

if(isFront) % TA, IP
    lm  = lmBar+ma*qD;
    dLm = ma*dqD;
else
    lm  = lmBar-ma*qD;
    dLm = -ma*dqD;
end
if(act<0)
    act=0;
end

xi  = lm/lmBar;
fmPE=0.0159*fmBar*exp(5.9*lmBar*(xi-1)-1); % Parallel elastic elements
k = 0.32+0.71*exp(-1.112*(xi-1))*sin(3.722*(xi-0.656));
eta = dLm/dLmBar;
h = 1+tanh(3.0*eta);
fm= fmBar*k*h*act+fmPE; 

if(fm<0)
    fm = 0;
elseif(fm>fmBar)
    fm = fmBar;
end
end

%% Muscle force by Stiffness Control
function fmStiff = calcFmStiff(lm, dLm, muscleName, isAnkle)
global bodyParam ctrlParam

lmBar  = bodyParam.([muscleName '_lmBar']);
ma     = bodyParam.([muscleName '_ma']);

if(isAnkle)
    kP = ((1/ma)^2)*ctrlParam.kPAStiff;
    kD = ((1/ma)^2)*ctrlParam.kDAStiff;
else
    kP = ((1/ma)^2)*ctrlParam.kPHStiff;
    kD = ((1/ma)^2)*ctrlParam.kDHStiff;
end

fm      = -kP*(lm-lmBar) -kD*dLm;
fmStiff = max([-fm,0]);
end
