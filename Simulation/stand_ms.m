function stand_ms()
clear all;
close all;

global bodyParam simParam ctrlParam


simParam.useReflection = 0; % add/not add impulsive GC activity after floor tilt onset 
if(~simParam.useReflection)
    outDir = 'simResults';
else
    outDir = 'simResultsReflection';
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
simParam.outDir = outDir;

%- Body Parameters (mainly from SuzukiJTheorBiol2012) -%
height = 1.7; % Body Height: 170 (cm)
mass   = 60;  % Body Mass 60 (kg)
% Lower Link
m1 = mass*0.35;    % (kg)   Mass of Lower link
L1 = height*0.51;  % (m)    Length of Lower link
h1 = height*0.255; % (m)    Length from ankle to Lower link COM
J1 = m1*(h1^2);    % (kgm2) Moment of inertia
% Upper Link
m2 = mass*0.62;    % (kg)   Mass of Upper link
L2 = height*0.45;  % (m)    Length of Upper link
h2 = height*0.225; % (m)    Length from hip to Upper link COM
J2 = m2*(h2^2);    % (kgm2) Moment of inertia
%
bodyParam.g = 9.8;
bodyParam.m1 = m1;
bodyParam.h1 = h1;
bodyParam.L1 = L1;
bodyParam.J1 = J1;
%
bodyParam.m2 = m2;
bodyParam.h2 = h2;
bodyParam.L2 = L2;
bodyParam.J2 = J2;
m   = m1+m2;
h   = (m1*h1 +m2*(L1+h2))/m;
mgh = m*bodyParam.g*h;
%--%

%- Muscle Parameters -%
% TA
bodyParam.TA_fmBar  = 1140 ;  % (N)   Maximul isometric contraction force
bodyParam.TA_lmBar  = 0.10;  % (m)   Natural length
bodyParam.TA_dLmBar = 3.0 ;  % (m/s) Maximum contraction velocity
bodyParam.TA_ma     = 0.042; % (m)   Moment arm (m)
% GC
bodyParam.GC_fmBar  = 1605 ; % (N)
bodyParam.GC_lmBar  = 0.11;  % (m)
bodyParam.GC_dLmBar = 3.0;   % (m/s)
bodyParam.GC_ma     = 0.057;  
% IP
bodyParam.IP_fmBar  = 1395 ; % (N)
bodyParam.IP_lmBar  = 0.10;  % (m)
bodyParam.IP_dLmBar = 3.0;   % (m/s)
bodyParam.IP_ma     = 0.011; 
% GMax
bodyParam.GMax_fmBar  = 1300 ; % (N)
bodyParam.GMax_lmBar  = 0.14;  % (m)
bodyParam.GMax_dLmBar = 3.0;   % (m/s)
bodyParam.GMax_ma     = 0.047; 
%
bodyParam.tiltMax = 4; % Maximum Floor tilt angle (degree)
%-------%

%- Simulation Initial States -%
simParam.dq10 =  0.0;  % Ankle Angular Vel
simParam.dq20 =  0.0;  % Hip Angular Vel
simParam.q10  = -0.11; % Ankle Angle (rad) From Experiment
simParam.q20  =  0.18; % Hip   Angle (rad) From Experiment
simParam.q0   = [simParam.dq10; simParam.q10; simParam.dq20; simParam.q20];
% Simulation Condition
simParam.simSamp = 1/1000; % Sampling Freqency of Simulation: 1(kHz)
simParam.outSamp = 1/100;  % Output Frequency of Simulation: 100(Hz)
simParam.delay   = 0.15;   % Sensory Delay (s)
simParam.tEnd    = 20;     % Duration of Simulation (s)
%--%

%- Parameters for Passive Stiffness Control -%
ctrlParam.kPAStiff = 0.3*mgh;  % Propotional control gain (Kp) for Stiffness control around Ankle
ctrlParam.kPHStiff = 0.4*mgh;  % Propotional control gain (Kp) for Stiffness control around Hip
%
ctrlParam.kDAStiff = 4.0;  % Derivative control gain (Kd) for Stiffness control around Ankle
ctrlParam.kDHStiff = 10.0; % Derivative control gain (Kd) for Stiffness control around Hip
%--%

ctrlParam.sigma   = 0.33;  % funatoRSOS20162016
ctrlParam.ref = calcCOMX(simParam.q0, [0,0]); % Target COM X Position=61mm
disp(['Target COM X Position: ' num2str(ctrlParam.ref*1000, '%.0f') ' (mm)']);
ctrlParam.mpcSamp = 0.1;   % Sampling Freqency of MPC: 10(Hz)

% Name of the mex file for the model predictive controller
ctrlParam.mpcFuncName = 'mpcCtrlFunc'; 

%- Execute simulation (calcMatlab function) -%
simParam.simNum = 20; % Number of the simulation trials with different random number for noise
[t, outData] = calcMatlab(simParam, ctrlParam.ref); 
figure;
showResult(t, outData, t(end)); % Display results
%--%
end


%% Simulation using Matlab MPC Toolbox
function [t,outData] = calcMatlab(simParam, ref)
global last_mv 
global onlineData

%- Simulation Parameters -%
sampling = simParam.outSamp; % Sampling frequency for save & output: 10 ms
tEnd = simParam.tEnd;

%== Calculation of Differential equation ==%
for simIndx = 1:simParam.simNum
    simParam.simIndx = simIndx;
    %- Setting random numbers for noise -%
    rng("shuffle");
    simParam.seed = randi(100000000); % set seed for random number as random (this "seed" is saved for future replay)
    rng(simParam.seed); % set seed of random number for noise
    %--%

    %- Simulation Condition-%
    evalFuncU = [0.01, 0.01, 0.01]; % weights for control inputs in evaluation function
    refMod    = ref;                % fix COM reference
    % Save Simulation Condition
    simInfo = array2table([refMod, simParam.seed, evalFuncU],'VariableNames',{'Ref','Seed', 'EvalFuncU(1)', 'EvalFuncU(2)', 'EvalFuncU(3)'});
    writetable(simInfo, [simParam.outDir '/simInfo_' num2str(simParam.simIndx, '%03d') '.csv']);
    %--%
    % Initialize Model
    last_mv = rand([3,1])/2; % set initial inputs (random value between 0-0.5)
    q0 =simParam.q0;         % set initial state
    % Build MEX file for MPC Controller
    nlobj = setMPCCtrl(q0, evalFuncU);
    [~, onlineData] = getCodeGenerationData(nlobj, [simParam.q0; 0.0; 0.0], last_mv); % Initialize
    onlineData.ref = refMod; % COM reference
    %--%
    
    %- Calculation: differential equation is "diffEqLink2_ms.m" -%
    disp(['----- Start Calculation, Trial Number: ' num2str(simParam.simIndx) '/' num2str(simParam.simNum) ' -----'])
    [t,outData] = odeEuler('diffEqLink2_ms', 0.0:sampling:tEnd, q0, 0.0:simParam.simSamp:tEnd, sampling/simParam.simSamp,simParam.delay, simParam.simIndx); % Euler
    %--%
end
%==%
end

%% Euler法（Euler-Maruyama Method）
function [t,outData] = odeEuler(dFuncName, t, x0, tCalc, simRate, delay, simIndx)
global  simOutML simParam

dFunc = str2func(dFuncName);
dt = tCalc(2)-tCalc(1);        % Simulation Interval: 1 ms
lagTimeIndx = round(delay/dt); % Number of lag with delay: 0.15/0.001=150 sample

outData.phi     = zeros(length(t),2);
outData.tauFB   = zeros(length(t),1);
outData.actMV   = zeros(length(t),3);
outData.muscleU = zeros(length(t),2);
outData.uStiff  = zeros(length(t),2);
outData.U       = zeros(length(t),2);
outData.x       = zeros(length(t),length(x0));
outData.comX    = zeros(length(t),1);
outData.headZ   = zeros(length(t),1);
%
outData.noise   = zeros(length(t),2);
outData.mEqDDQ  = zeros(length(t),2);

xCalc = zeros(length(x0), length(tCalc));
xCalc(:,1)=x0; % Initial value

num=1;
tDispNext = 0;
for i=2:length(tCalc)
    if(tCalc(i)>tDispNext)
        disp(['time: ' num2str(tCalc(i), '%.1f')]);
        tDispNext = tCalc(i)+0.1; % display status every 0.1 (s)
    end    
    if(i > lagTimeIndx+1)
        z=xCalc(:,i-1-lagTimeIndx); 
    else
        z=x0;
    end
    xCalc(:,i) = xCalc(:,i-1)+dFunc(tCalc(i),xCalc(:,i-1), z)*dt;

    if(i==2||rem(i,simRate)==0)
        outData.phi(num,:)     = simOutML.phiList;
        outData.actMV(num,:)   = simOutML.mv;
        outData.muscleU(num,:) = simOutML.muscleU;
        outData.uStiff(num,:)  = simOutML.uStiff;
        outData.U(num,:)       = simOutML.U;
        outData.x(num,:)       = xCalc(:,i)';
        %
        outData.comX(num, 1)   = calcCOMX(xCalc(:,i), simOutML.phiList); 
        outData.headZ(num,1)   = calcHeadZ(xCalc(:,i), simOutML.phiList);
        outData.noise(num,:)   = simOutML.noise;
        outData.mEqDDQ(num,:)  = simOutML.mEqDDQ;
        %
        dispProgress = 1;
        if(dispProgress)
            if(i==2)
                figure;
            else
                showResult(t, outData, t(num));
                drawnow
            end
        end
        %- Output simulation Results: "simResultData_{trial number}.csv" -%
        saveResult =1;
        if(saveResult)
            resultFilename = [simParam.outDir '/simResultData_' num2str(simIndx, '%03d') '.csv'];
            saveData(resultFilename, t, outData, num);         
        end
        %--%
        num=num+1;
    end
end
end

%% Setting Model Predictive Controller
function nlobj = setMPCCtrl(q0, evalFuncU)
global ctrlParam last_mv
global currentTime 

% Structure of the Controller
nx = 6; % Number of states: dqD1,dqD2,qD1,qD2, t, phi
ny = 1; % Number of outputs: comX
nu = 3; % Number of inputs: mv1(TA),mv2(GC), mv3(IP/GM)
nlobj = nlmpc(nx, ny, nu);
% Sampling rate of Prediction and Controller 
nlobj.Ts = ctrlParam.mpcSamp;
% Prediction horizon
nlobj.PredictionHorizon = 30; % good
% Control horizon
nlobj.ControlHorizon = 2;

% Nonlinear Plant Model
nlobj.Model.StateFcn = "odeEulerModel"; % Internal model: The contents are described in "odeEulerModel.m".
nlobj.Model.IsContinuousTime = false;   % (Required to use explicitly specify ODE instead of using default ODE)
nlobj.Model.NumberOfParameters = 0;
% Solver setting
nlobj.Optimization.SolverOptions.Display = 'none';

% Output function (output: comX): The contents are described in "standOutFun_ms.m"
nlobj.Model.OutputFcn = 'standOutFunc_ms'; % 

% Weights for the evaluation fucntion
nlobj.Weights.OutputVariables = 1.0;                % Weight for Output (comX)
nlobj.Weights.ManipulatedVariables = evalFuncU;     % Weight for Inputs
nlobj.Weights.ManipulatedVariablesRate = [1, 1, 1]; % Wdight for Input Changes

% Constraints: Range of the muscle activities
nlobj.MV(1).Min =  0; % TA 
nlobj.MV(2).Min =  0; % GC
nlobj.MV(3).Min = -1; % IP/GM
nlobj.MV(1).Max =  1;
nlobj.MV(2).Max =  1;
nlobj.MV(3).Max =  1;
% Set the above constraints as hard constraints
nlobj.MV(1).MinECR = 0; % hard constraints
nlobj.MV(1).MaxECR = 0; % hard constraints
nlobj.MV(2).MinECR = 0; % hard constraints
nlobj.MV(2).MaxECR = 0; % hard constraints
nlobj.MV(3).MinECR = 0; % hard constraints
nlobj.MV(3).MaxECR = 0; % hard constraints
% Constraints: range of the output (comX)
nlobj.OV(1).Min = -0.10; % (m) maximum posterior com position: 100mm
nlobj.OV(1).Max =  0.20; % (m) maximum anterior com position:  200mm

%- Build MEX file -%
t0 = 0.0; % Initial time
currentTime = t0;
x0 = [q0; t0; 0]; % Set initial time and states
[coreData,onlineData] = getCodeGenerationData(nlobj,x0,last_mv);
onlineData.ref = ctrlParam.ref; % Set reference (comX)
buildMEX(nlobj, ctrlParam.mpcFuncName ,coreData, onlineData); % Build MEX function
%--%

end

%% Display Results
function showResult(t, outData, tEnd)
q = outData.x;

subplot(4,4,1)
plot(t,q(:,2));
ylabel('theta ankle (rad)');
xlim([t(1),tEnd]);
%
subplot(4,4,5)
plot(t,q(:,4));
ylabel('theta hip (rad)');
xlim([t(1),tEnd]);
%
subplot(4,4,9)
plot(t,outData.comX*1000); 
ylabel('COM X (mm)');
xlim([t(1),tEnd]);
%
subplot(4,4,13)
plot(t,outData.headZ);
ylabel('Head Z (m)');
xlim([t(1),tEnd]);
%
subplot(4,4,2)
plot(t,outData.actMV(:,1));
ylabel('Activation TA');
xlim([t(1),tEnd]);
%
subplot(4,4,6)
plot(t,outData.actMV(:,2));
ylabel('Activation GAS');
xlim([t(1),tEnd]);
%
subplot(4,4,10)
plot(t,max(outData.actMV(:,3),0));
ylabel('Activation IP');
xlim([t(1),tEnd]);
%
subplot(4,4,14)
plot(t,max(-outData.actMV(:,3),0));
ylabel('Activation GM');
xlim([t(1),tEnd]);
%
subplot(4,4,3)
plot(t,outData.muscleU(:,1));
ylabel('Ankle Muscle Torque');
xlim([t(1),tEnd]);
%
subplot(4,4,7)
plot(t,outData.muscleU(:,2));
ylabel('Hip Muscle Torque');
xlim([t(1),tEnd]);
%
subplot(4,4,4)
plot(t,outData.phi(:,1));
ylabel('floorAng (rad)');
xlim([t(1),tEnd]);
%
subplot(4,4,8)
plot(t,outData.phi(:,2));
ylabel('floorVel (rad/s)');
xlim([t(1),tEnd]);
%
subplot(4,4,11)
plot(t,outData.uStiff(:,1));
ylabel('Muscle Stiff Ankle');
xlim([t(1),tEnd]);
%
subplot(4,4,15)
plot(t,outData.uStiff(:,2));
ylabel('Muscle Stiff Hip');
xlim([t(1),tEnd]);
%
subplot(4,4,12)
plot(t,outData.U(:,1));
ylabel('Total Ankle Torque');
xlim([t(1),tEnd]);
%
subplot(4,4,16)
plot(t,outData.U(:,2));
ylabel('Total Hip Torque');
xlim([t(1),tEnd]);
end

%% Calculate Horizontal COM Position (comX)
function comX = calcCOMX(x, phiList)
global bodyParam

phi  = phiList(1);
m    = bodyParam.m1+bodyParam.m2;
comLegX   = -bodyParam.h1*sin(x(2)+phi);
comTrunkX = -bodyParam.L1*sin(x(2)+phi)-bodyParam.h2*sin(x(4)+x(2)+phi);
comX = (bodyParam.m1*comLegX + bodyParam.m2*comTrunkX)/m;

end

%% Calculate Head Height
function headZ = calcHeadZ(x, phiList)
global bodyParam

phi  = phiList(1);
headZ = bodyParam.L1*cos(x(2)+phi)+bodyParam.L2*cos(x(4)+x(2)+phi);
end

%% Save Results
function saveData(resultFilename, t, outData, num)
if(num==1)
    text ={'time', 'dqAnkle', 'qAnkle', 'dqHip', 'qHip', 'ActivationTA','ActivationGAS', ...
        'ActivationIP', 'ActivationGM', 'AnkleTorque', 'HipTorque', 'FloorAngle', 'FloorVel', ...
        'StiffnessAnkle', 'StiffnessHip', 'COM X', 'Head Z', 'Total Ankle Torque', 'Total Hip Torque'};
    fid = fopen(resultFilename, 'w','n','Shift_JIS');
    for tI = 1:(length(text)-1)
        fprintf(fid, '%s, ', text{tI});
    end
    fprintf(fid, '%s\n', text{end});
    fclose(fid);
end
fid = fopen(resultFilename, 'A');
fprintf(fid, '%s, ', num2str(t(num)));
fprintf(fid, '%s, ', num2str(outData.x(num,1)));
fprintf(fid, '%s, ', num2str(outData.x(num,2)));
fprintf(fid, '%s, ', num2str(outData.x(num,3)));
fprintf(fid, '%s, ', num2str(outData.x(num,4)));
fprintf(fid, '%s, ', num2str(outData.actMV(num,1)));
fprintf(fid, '%s, ', num2str(outData.actMV(num,2)));
fprintf(fid, '%s, ', num2str(outData.actMV(num,3)));
fprintf(fid, '%s, ', num2str(-outData.actMV(num,3)));
fprintf(fid, '%s, ', num2str(outData.muscleU(num,1)));
fprintf(fid, '%s, ', num2str(outData.muscleU(num,2)));
fprintf(fid, '%s, ', num2str(outData.phi(num,1)));
fprintf(fid, '%s, ', num2str(outData.phi(num,2)));
fprintf(fid, '%s, ', num2str(outData.uStiff(num,1)));
fprintf(fid, '%s, ', num2str(outData.uStiff(num,2)));
fprintf(fid, '%s, ', num2str(outData.comX(num,1)));
fprintf(fid, '%s, ', num2str(outData.headZ(num,1)));
fprintf(fid, '%s, ', num2str(outData.U(num,1)));
fprintf(fid, '%s\n', num2str(outData.U(num,2)));
fclose(fid);
end

