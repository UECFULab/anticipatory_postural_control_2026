function showResultFromData()
clear all;
close all;

showReflectionResult = 0;
saveFigs = 0;

if(~showReflectionResult)
    dataDir = 'simResults';
    simNum = 20;
    figDir ='./figs/Fig6_sim/';
else
    dataDir = 'simResultsReflection';
    simNum = 20;
    figDir ='./figs/Fig9_simReflection/';
end
if(saveFigs)
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
end
%
[timeAll, COMAll, taAll, gcAll, ankAngAll, hipAngAll] = loadFile(dataDir, simNum);
time = timeAll(:,1)-15; % Set floor tilt start time (15s) as 0s

if(showReflectionResult)
    sFreq  = 100;
    lpDim  = 4;
    lpFreq = 10;
    [b,a] = butter(lpDim, (lpFreq*2)/sFreq, 'low');
    for j=1:size(gcAll,2)
        taAll(:,j) = filtfilt(b,a,taAll(:,j));
        gcAll(:,j) = filtfilt(b,a,gcAll(:,j));
    end
end

%- Time series of COM X -%
figure;
subplot(2,1,1)
plot(time, mean(COMAll(:,1:end)*1000,2)); hold on;
plot(time, mean(COMAll(:,1:end)*1000,2)+std(COMAll(:,1:end)*1000,0,2), 'r-');
plot(time, mean(COMAll(:,1:end)*1000,2)-std(COMAll(:,1:end)*1000,0,2), 'r-'); 
ylim([0, 100])
yticks(0:50:100)
xline(0,'r');
xline(1,'r');
xline(-2,'b');
hold off;
ylabel('COM (mm)')
xlim([-5,5])
if(saveFigs)
    saveas(gcf, [figDir 'FigA_simCOM'],'emf')
end
%--%

%- Time series of Hip and Ankle -%
figure;
subplot(2,1,1)
plot(time, mean(hipAngAll(:,1:end),2)); hold on;
plot(time, mean(hipAngAll(:,1:end),2)+std(hipAngAll(:,1:end),0,2), 'r-');
plot(time, mean(hipAngAll(:,1:end),2)-std(hipAngAll(:,1:end),0,2), 'r-'); 
if(~showReflectionResult)
    ylim([-0.1,0.2])
    yticks(-0.1:0.1:0.2)
else
    ylim([-0.1,0.3])
    yticks(-0.2:0.2:0.4)
end
xlim([-5,5])
xline(0,'r');
xline(1,'r');
xline(-2,'b');
hold off;
ylabel('Hip (rad)')
%
subplot(2,1,2)
plot(time, mean(ankAngAll(:,1:end),2)); hold on;
plot(time, mean(ankAngAll(:,1:end),2)+std(ankAngAll(:,1:end),0,2), 'r-');
plot(time, mean(ankAngAll(:,1:end),2)-std(ankAngAll(:,1:end),0,2), 'r-'); 
ylim([-0.3, 0.0])
yticks(-0.3:0.1:0.0)
xlim([-5,5])
xline(0,'r');
xline(1,'r');
xline(-2,'b');
hold off;
ylabel('Ankle (rad)')
if(saveFigs)
    saveas(gcf, [figDir 'FigB_simAnkHip'],'emf')
end
%--%

%- Time series of TA and GC -% 
figure;
subplot(2,1,1)
plot(time, mean(taAll(:,1:end),2)*100); hold on;
plot(time, mean(taAll(:,1:end),2)*100+std(taAll(:,1:end)*100,0,2), 'r-');
plot(time, mean(taAll(:,1:end),2)*100-std(taAll(:,1:end)*100,0,2), 'r-'); 
xline(0,'r');
xline(1,'r');
xline(-2,'b');
hold off;
ylabel('TA')
ylim([0,100])
yticks([0,50,100])
ylim([0,100])
yticks([0,50,100])
xlim([-5,5])
subplot(2,1,2)
plot(time, mean(gcAll(:,1:end),2)*100); hold on;
plot(time, mean(gcAll(:,1:end),2)*100+std(gcAll(:,1:end)*100,0,2), 'r-');
plot(time, mean(gcAll(:,1:end),2)*100-std(gcAll(:,1:end)*100,0,2), 'r-'); 
xline(0,'r');
xline(1,'r');
xline(-2,'b');
hold off;
ylabel('GC')
ylim([0,100])
yticks([0,50,100])
xlim([-5,5])
if(saveFigs)
    saveas(gcf, [figDir 'FigE_simTAGC'],'emf')
end
%--%

soundStartIndx = 13*100+1;
floorStartIndx = 15*100+1;


%== Fig CS: COM/Hip/Ankle Bar plot ==%
%- Bar plot: COM/Hip/Ankle values at CS and FS -%
figure;
subplot(3,1,1)
drawErrorbarMSD([mean(COMAll(soundStartIndx,:)*1000), mean(COMAll(floorStartIndx,:)*1000)], [std(COMAll(soundStartIndx,:)*1000), std(COMAll(floorStartIndx,:)*1000)]);
ylim([0, 120])
yticks(0:50:100)
ylabel('COM');
xticklabels({'CS','FS'})
%
subplot(3,1,2)
drawErrorbarMSD([mean(hipAngAll(soundStartIndx,:)), mean(hipAngAll(floorStartIndx,:))], [std(hipAngAll(soundStartIndx,:)), std(hipAngAll(floorStartIndx,:))]);
ylim([0, 0.25])
yticks(0:0.1:0.2)
ylabel('Hip');
xticklabels({'CS','FS'})
%
subplot(3,1,3)
drawErrorbarMSD([mean(ankAngAll(soundStartIndx,:)), mean(ankAngAll(floorStartIndx,:))], [std(ankAngAll(soundStartIndx,:)), std(ankAngAll(floorStartIndx,:))]);
ylim([-0.2, 0])
yticks(-0.2:0.1:0)
ylabel('Ankle');
xticklabels({'CS','FS'})
if(saveFigs)
    saveas(gcf, [figDir 'FigCDLeft_simCSFSCOMHipAnk'],'emf')
end
% t-test between CS and FS
% COM
[~,pCOM, ~, statesCOM] = ttest2(COMAll(soundStartIndx,:), COMAll(floorStartIndx,:));
disp(['COM diff: P=' num2str(pCOM) ', t-value=' num2str(statesCOM.tstat) ', dof=' num2str(statesCOM.df)]);
% Ankle
[~,pHip, ~, statesHip] = ttest2(hipAngAll(soundStartIndx,:), hipAngAll(floorStartIndx,:));
disp(['Hip diff: P=' num2str(pHip) ', t-value=' num2str(statesHip.tstat) ', dof=' num2str(statesHip.df)]);
% Hip
[~,pAnk, ~, statesAnk] = ttest2(ankAngAll(soundStartIndx,:), ankAngAll(floorStartIndx,:));
disp(['Ankle diff: P=' num2str(pAnk) ', t-value=' num2str(statesAnk.tstat) ', dof=' num2str(statesAnk.df)]);
%--%
%
%- Diff CS-FS -%
COMDiffAll   = COMAll(floorStartIndx,:)-COMAll(soundStartIndx,:);
hipDiffAll   = hipAngAll(floorStartIndx,:)-hipAngAll(soundStartIndx,:);
ankleDiffAll = ankAngAll(floorStartIndx,:)-ankAngAll(soundStartIndx,:);
%
figure;
subplot(3,1,1)
drawErrorbarMSD([mean(COMDiffAll*1000)], [std(COMDiffAll*1000)]);
ylim([0, 30])
yticks(0:15:30)
ylabel('COM');
%
subplot(3,1,2)
drawErrorbarMSD([mean(hipDiffAll)], [std(hipDiffAll)]);
ylim([-0.04, 0])
yticks(-0.04:0.02:0)
ylabel('Hip');
%
subplot(3,1,3)
drawErrorbarMSD([mean(ankleDiffAll)], [std(ankleDiffAll)]);
ylim([-0.04, 0])
yticks(-0.04:0.02:0)
ylabel('Ankle');
if(saveFigs)
    saveas(gcf, [figDir 'FigCDRight_simDiffCOMHipAnk'],'emf')
end
%--%
%== ==%

%== Fig F: EMG Bar plot ==%
%- Bar plot: TA/GC values at CS and FS -%
figure;
subplot(2,1,1)
drawErrorbarMSD([mean(taAll(soundStartIndx,:)*100), mean(taAll(floorStartIndx,:)*100)], [std(taAll(soundStartIndx,:)*100), std(taAll(floorStartIndx,:)*100)]);
ylim([0, 50])
yticks(0:25:50)
ylabel('TA');
xticklabels({'CS','FS'})
%
subplot(2,1,2)
drawErrorbarMSD([mean(gcAll(soundStartIndx,:)*100), mean(gcAll(floorStartIndx,:)*100)], [std(gcAll(soundStartIndx,:)*100), std(gcAll(floorStartIndx,:)*100)]);
ylim([0, 50])
yticks(0:25:50)
ylabel('GC');
xticklabels({'CS','FS'})
if(saveFigs)
    saveas(gcf, [figDir 'FigFLeft_simCSFSEMG'],'emf')
end
% t-test between CS and FS
% TA
[~,pTA, ~, statesTA] = ttest2(taAll(soundStartIndx,:), taAll(floorStartIndx,:));
disp(['TA diff: P=' num2str(pTA) ', t-value=' num2str(statesTA.tstat) ', dof=' num2str(statesTA.df)]);
% GC
[~,pGC, ~, statesGC] = ttest2(gcAll(soundStartIndx,:), gcAll(floorStartIndx,:));
disp(['GC diff: P=' num2str(pGC) ', t-value=' num2str(statesGC.tstat) ', dof=' num2str(statesGC.df)]);
%--%
%
%- Grad CS-FS -%
taCoeffGroup = [];
gcCoeffGroup = [];
for i=1:size(taAll,2)
    taCoeffGroup(i) = linearRegressionMatlab(time(soundStartIndx:floorStartIndx)', taAll(soundStartIndx:floorStartIndx,i)'*100);
    gcCoeffGroup(i) = linearRegressionMatlab(time(soundStartIndx:floorStartIndx)', gcAll(soundStartIndx:floorStartIndx,i)'*100);
end
%
figure;
subplot(2,1,1)
drawErrorbarMSD([mean(taCoeffGroup)], [std(taCoeffGroup)]);
% ylim([-1, 6])
% yticks(0:2:6)
ylim([-1, 4])
yticks(0:2:4)
ylabel('TA');
%
subplot(2,1,2)
drawErrorbarMSD([mean(gcCoeffGroup)], [std(gcCoeffGroup)]);
% ylim([-1, 6])
% yticks(0:2:6)
ylim([-1, 4])
yticks(0:2:4)
ylabel('GC');
if(saveFigs)
    saveas(gcf, [figDir 'FigFRight_simEMGCoeff'],'emf')
end
%== ==%

if(~showReflectionResult)
    %== Fig G: Relationship between COM Diff and GC coefficient ==%
    [cVal, r2Val, pVal, ~, bVal, ~] = linearRegressionMatlab(COMDiffAll*1000, gcCoeffGroup);
    disp(['R2 = ' num2str(r2Val) ', P=' num2str(pVal) ', c=' num2str(cVal)]);
    figure;
    hold on;
    plot(COMDiffAll*1000, gcCoeffGroup, '.')
    % xlim([15,25])
    % xticks(15:5:25)
    % ylim([-3, 6]);
    % yticks(-3:3:6)
    xlim([15,30])
    xticks(15:5:30)
    ylim([-4, 6]);
    yticks(-3:3:6)
    xl = xlim;
    xVal = xl(1):0.05:xl(2);
    plot(xVal, cVal*xVal+bVal, '-');
    hold off;
    title('COM-GC(coeff)')
    if(saveFigs)
        saveas(gcf, [figDir 'FigG_simCOMGC'],'emf')
    end
    %== ==%

    %== Fig H: COM phase portrait ==%
    % calculate COM velocity
    mFreq    = 100; % OutputFrequency: 100 (Hz)
    interval = 1/mFreq;
    lpDim  = 4;
    lpFreq = 5;
    [b,a] = butter(lpDim, (lpFreq*2)/ mFreq, 'low');
    dCOMAll = zeros(size(COMAll,1)-4, size(COMAll,2));
    for i=1:size(COMAll,2)
        dCOMAll(:,i) = diff5p(filtfilt(b,a,COMAll(:,i)),  interval);
    end
    %
    COMpredAll  = COMAll(soundStartIndx:floorStartIndx,:);
    dCOMpredAll = dCOMAll(soundStartIndx+2:floorStartIndx+2,:);
    GCpredAll   = gcAll(soundStartIndx:floorStartIndx,:)*100;
    %
    % Down sample (for reducing the output file size)
    COMpredDS  = downsample(COMpredAll, 2);
    dCOMpredDS = downsample(dCOMpredAll, 2);
    GCpredDS   = downsample(GCpredAll, 2);
    %
    % Discretize the GCpredDS value into 20 levels
    discNum = 20;
    edges = linspace(0, max(GCpredDS(:)), discNum+1);
    GCDS_discrete = discretize(GCpredDS, edges);
    %
    figure;
    cm = colormap(nebula(discNum));
    hold on;
    for triIndx = 1:size(COMpredDS,2)
        for i=1:(size(COMpredDS,1)-1)
            colorVal = cm(GCDS_discrete(i, triIndx),:);
            plot(COMpredDS(i:i+1, triIndx)*1000, dCOMpredDS(i:i+1, triIndx)*1000, '-', 'Color',colorVal, 'LineWidth',6);
        end
    end
    hold off;
    axis square;
    clim([0 max(GCpredDS(:))])
    xlabel('COM (mm)')
    ylabel('COM Velocity(mm/s)')
    ylim([-20,60])
    set(gca, 'ytick', -20:20:60);
    %
    colorbar;
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 42)
    if(saveFigs)
        exportgraphics(gcf, [figDir 'FigH_COM-dCOM_sim.emf'], 'ContentType', 'vector');
    end
    %== ==%
end
end

%%
function drawErrorbarMSD(dataM, dataSD)
bar(dataM);     
hold on;
eh0=errorbar(1:length(dataM), dataM, dataSD);
hold off;
xlim([0.5 length(dataM)+0.5]);
set(gca, 'xtick', 1:1:length(dataM));
set(eh0, 'LineStyle', 'none', 'LineWidth', 2);
end

%%
function [timeAll, COMAll, taAll, gcAll, ankAngAll, hipAngAll] = loadFile(dataDir, simNum)

for simIndx = 1:simNum
    resultFilename = [dataDir '/simResultData_' num2str(simIndx, '%03d') '.csv'];
    M=readmatrix(resultFilename);
    resultData.time   = M(:,1);
    resultData.x      = M(:,2:5);
    resultData.comX   = M(:,16);
    resultData.mvTA = M(:,6);  % TA
    resultData.mvGAS = M(:,7); % GC 
    if(~exist('timeAll'))
        timeAll= zeros(length(resultData.time), simNum);
        COMAll = zeros(length(resultData.time), simNum);
        taAll  = zeros(length(resultData.time), simNum);
        gcAll  = zeros(length(resultData.time), simNum);
        ankAngAll  = zeros(length(resultData.time), simNum);
        hipAngAll  = zeros(length(resultData.time), simNum);
    end
    timeAll(:,simIndx)= resultData.time;
    COMAll(:,simIndx) = resultData.comX;
    taAll(:,simIndx)  = resultData.mvTA;
    gcAll(:,simIndx)  = resultData.mvGAS;
    ankAngAll(:,simIndx)  = resultData.x(:,2);
    hipAngAll(:,simIndx)  = resultData.x(:,4);
end
end

%% Linear Regression
function [cVal, r2Val, pVal, rmseVal, bVal, lm] = linearRegressionMatlab(xSeries, ySeries)

tbl = table(xSeries',ySeries','VariableNames',{'xSeries','ySeries'});
lm = fitlm(tbl, 'ySeries~xSeries','RobustOpts','on'); 
coeff = lm.Coefficients;
cVal   = coeff{'xSeries', 'Estimate'};
r2Val  = lm.Rsquared.Ordinary;
pVal   = coeff{'xSeries', 'pValue'};
rmseVal= lm.RMSE;
bVal   = coeff{'(Intercept)', 'Estimate'};
end

%% 5-point finite difference
function df = diff5p(dat,interval)

n=1;
df=zeros(1, length(dat)-4);
for i=3:(length(dat)-2)
    df(n)=(dat(i-2)-8*dat(i-1)+8*dat(i+1)-dat(i+2))/(12*interval);
    n=n+1;
end
end