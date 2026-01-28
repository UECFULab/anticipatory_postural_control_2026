function analysisEMG()
clear all
close all

matDir = 'data/';
figDir = 'figs/';
saveFigs = 1; 
if(saveFigs)
    if ~exist(figDir, 'dir')
        mkdir(figDir);
    end
end

extractDur = 5; % Duration (s) for analysis with Floor tilt start as 0 (s)
sList = 1:10;
% sList = 1:2;

cTAnonSList  = cell(1,length(sList));
cTAconSList  = cell(1,length(sList));
cGCnonSList = cell(1,length(sList));
cGCconSList = cell(1,length(sList));
%
cTAnonSListM  = zeros(1, length(sList));
cTAconSListM  = zeros(1, length(sList));
cGCnonSListM = zeros(1, length(sList));
cGCconSListM = zeros(1, length(sList));
%
pTAnonSListM  = zeros(1, length(sList));
pTAconSListM  = zeros(1, length(sList));
pGCnonSListM = zeros(1, length(sList));
pGCconSListM = zeros(1, length(sList));
for sIndx = sList
    Name = ['subject' num2str(sIndx)]; % Subject Number
    if(strcmp(Name, 'subject1'))
        SkipData = 2;
    elseif(strcmp(Name, 'subject2'))
        SkipData = [];
    elseif(strcmp(Name, 'subject3'))
        SkipData = [1, 57]; % Motion: 1, 57; encoder: 42, 56
    elseif(strcmp(Name, 'subject4'))
        SkipData = [];
    elseif(strcmp(Name, 'subject5'))
        SkipData = 10; 
    elseif(strcmp(Name, 'subject6'))
        SkipData = [16,20];
    elseif(strcmp(Name, 'subject7'))
        SkipData = 1; 
    elseif(strcmp(Name, 'subject8'))
        SkipData = [];
    elseif(strcmp(Name, 'subject9'))
        SkipData = 44;
    elseif(strcmp(Name, 'subject10'))
        SkipData = 1; 
    else
        SkipData = [];
    end
    disp(['-- ', Name, '--'])

    %-- Read emg data --% 
    emgFileName = [matDir, 'EMGData_',Name,'.xlsx'];
    RTA = readmatrix(emgFileName, 'sheet', 'RTA');
    LTA = readmatrix(emgFileName, 'sheet', 'LTA');
    RGC = readmatrix(emgFileName, 'sheet', 'RGC');
    LGC = readmatrix(emgFileName, 'sheet', 'LGC');
    %--%

    %- Read motion data -%
    motionFileName = [matDir, 'MotionData_',Name,'.xlsx'];
    tiltEndIndx = readmatrix(motionFileName, 'sheet', 'tiltEndIndx'); % measured freqeuency: 300Hz
    Soundstart  = readmatrix(motionFileName, 'sheet', 'Soundstart');  % measured freqeuency: 1000Hz
    Soundend    = readmatrix(motionFileName, 'sheet', 'Soundend');    % measured freqeuency: 1000Hz
    %--%

    nonIndx = 11:30; 
    conIndx = 41:60; 
    nonIndx = nonIndx(~ismember(nonIndx, SkipData)); % Remove failed trial (SkipData)
    conIndx = conIndx(~ismember(conIndx, SkipData)); % Remove failed trial (SkipData)
    % TA
    RTANonS = RTA(:,nonIndx);
    LTANonS = LTA(:,nonIndx);
    RTAConS = RTA(:,conIndx);
    LTAConS = LTA(:,conIndx);
    TAnonS  = (RTANonS+LTANonS)/2;
    TAconS  = (RTAConS+LTAConS)/2;
    % GC
    RGCNonS = RGC(:,nonIndx);
    LGCNonS = LGC(:,nonIndx);
    RGCConS = RGC(:,conIndx);
    LGCConS = LGC(:,conIndx);
    GCnonS  = (RGCNonS+LGCNonS)/2;
    GCconS  = (RGCConS+LGCConS)/2;
    %
    EMGTime = linspace(-extractDur, extractDur, 2*extractDur*300);
    Soundstartmean = mean(Soundstart(conIndx),2)/1000;
    Soundendmean   = mean(Soundend(conIndx),2)/1000;
    SoundstartSD   = std(Soundstart(conIndx),[], 2)/1000;
    SoundendSD     = std(Soundend(conIndx),[], 2)/1000;
    %
    tiltEndMean    = mean(tiltEndIndx(conIndx),  2)/300;
    disp(['Cue start time (mean) : ' num2str(Soundstartmean, '%.2f')  ') (s), Cue end time (mean) : ' num2str(Soundendmean) ' (s), Tilt end time (mean) : ' num2str(tiltEndMean) ' (s)']);

    %- Plot EMG Figures for Each Subject -%
    %- No Cue -%
    figure;
    subplot(2,1,1); plotEMGRLFigs(EMGTime, TAnonS, [], [],  [], [], tiltEndMean, 'TA RL No Sound');
    subplot(2,1,2); plotEMGRLFigs(EMGTime, GCnonS, [], [],  [], [], tiltEndMean, 'GC RL No Sound');
    if(saveFigs)
        saveas(gcf, [figDir 'EMGNonS_' Name],'emf')
    end
    %--%
    %- With Cue -%
    figure;
    subplot(2,1,1); plotEMGRLFigs(EMGTime, TAconS, Soundstartmean, Soundendmean,  SoundstartSD, SoundendSD, tiltEndMean, 'TA RL With Sound');
    subplot(2,1,2); plotEMGRLFigs(EMGTime, GCconS, Soundstartmean, Soundendmean,  SoundstartSD, SoundendSD, tiltEndMean, 'GC RL With Sound');
    if(saveFigs)
        saveas(gcf, [figDir 'EMGConS_' Name],'emf')
    end
    %--%

    [~,soundStartIndx] = min(abs(EMGTime - Soundstartmean));
    [~,floorStartIndx] = min(abs(EMGTime - 0));
    EMGTimePreparation = EMGTime(soundStartIndx:floorStartIndx);
    TAnonSpreparation  = TAnonS(soundStartIndx:floorStartIndx,:);
    TAconSpreparation  = TAconS(soundStartIndx:floorStartIndx,:);
    GCnonSpreparation = GCnonS(soundStartIndx:floorStartIndx,:);
    GCconSpreparation = GCconS(soundStartIndx:floorStartIndx,:);
    %
    cTAnonS = zeros(1,size(GCconSpreparation,2));
    cTAconS = zeros(1,size(GCconSpreparation,2));
    cGCnonS = zeros(1,size(GCconSpreparation,2));
    cGCconS = zeros(1,size(GCconSpreparation,2));
    for i=1:size(GCnonSpreparation,2)
        cTAnonS(i) = linearRegressionMatlab(EMGTimePreparation, TAnonSpreparation(:,i)');
        cGCnonS(i) = linearRegressionMatlab(EMGTimePreparation, GCnonSpreparation(:,i)');
    end
    for i=1:size(GCconSpreparation,2)
        cTAconS(i) = linearRegressionMatlab(EMGTimePreparation, TAconSpreparation(:,i)');
        cGCconS(i) = linearRegressionMatlab(EMGTimePreparation, GCconSpreparation(:,i)');
    end
    cTAnonSList{sIndx}  = cTAnonS;
    cTAconSList{sIndx}  = cTAconS;
    cGCnonSList{sIndx} = cGCnonS;
    cGCconSList{sIndx} = cGCconS;
    %
    [cTAnonSListM(sIndx), ~, pTAnonSListM(sIndx)] = linearRegressionMatlab(EMGTimePreparation, mean(TAnonSpreparation, 2)');
    [cTAconSListM(sIndx), ~, pTAconSListM(sIndx)] = linearRegressionMatlab(EMGTimePreparation, mean(TAconSpreparation, 2)');
    [cGCnonSListM(sIndx), ~, pGCnonSListM(sIndx)] = linearRegressionMatlab(EMGTimePreparation, mean(GCnonSpreparation, 2)');
    [cGCconSListM(sIndx), ~, pGCconSListM(sIndx)] = linearRegressionMatlab(EMGTimePreparation, mean(GCconSpreparation, 2)');
end
%
figure;
subplot(2,1,1)
plotErrorbarGroup([cellfun(@mean, cTAnonSList)', cellfun(@mean, cTAconSList)'], [cellfun(@std, cTAnonSList)', cellfun(@std, cTAconSList)']);
ylim([-2 2])
set(gca, 'ytick', -2:1:2);
ylabel('TA')
subplot(2,1,2)
plotErrorbarGroup([cellfun(@mean, cGCnonSList)', cellfun(@mean, cGCconSList)'], [cellfun(@std, cGCnonSList)', cellfun(@std, cGCconSList)']);
ylabel('GC')
ylim([-10 12])
set(gca, 'ytick', -10:10:10);
if(saveFigs)
    saveas(gcf, [figDir 'cTAGC_ALL'],'emf')
end

%-- t-Test: difference in coefficient of linear regression between CS-FS --%
tResultTA = zeros(length(sList),3);
tResultGC = zeros(length(sList),3);
for sIndx = sList
    [~,tResultTA(sIndx,1),~,statesTA] = ttest2(cTAnonSList{sIndx}, cTAconSList{sIndx});
    [~,tResultGC(sIndx,1),~,statesGC] = ttest2(cGCnonSList{sIndx}, cGCconSList{sIndx});
    tResultTA(sIndx,2) = statesTA.tstat;
    tResultGC(sIndx,2) = statesGC.tstat;
    tResultTA(sIndx,3) = statesTA.df;
    tResultGC(sIndx,3) = statesGC.df;
end
tResultTATable = array2table(tResultTA,'VariableNames',{'p-value','t-value','df'});
tResultGCTable = array2table(tResultGC,'VariableNames',{'p-value','t-value','df'});
writetable(tResultTATable, 'tResultTA.csv');
writetable(tResultGCTable, 'tResultGC.csv');
%--%
end

%% Plot EMG Figures
function plotEMGRLFigs(EMGTime, emgData, Soundstartmean, Soundendmean,  SoundstartSD, SoundendSD, tiltEndMean, emgName)

hold on
emgMean = mean(emgData,2);
emgSD   = std(emgData,0,2);

plot(EMGTime, emgMean, 'b-','LineWidth',2)
plot(EMGTime, emgMean+emgSD, 'c-')
plot(EMGTime, emgMean-emgSD, 'c-')
h_axes = gca;
h_axes.FontSize = 12;
xlim([-5 5]);
set(gca, 'xtick', -5:1:5);
ylim([0 40])
if(~isempty(Soundstartmean))
    xline(Soundstartmean,'r-');
    xline(Soundstartmean+SoundstartSD,'m-');
    xline(Soundstartmean-SoundstartSD,'m-');
    xline(Soundendmean,'r-');
    xline(Soundendmean+SoundendSD,'m-');
    xline(Soundendmean-SoundendSD,'m-');
end
xline(0,'k','LineWidth',2);
xline(tiltEndMean,'k','LineWidth',2);
fontsize(10,"points")
xlabel('Time(s)')
ylabel('MVC (%)')
title(emgName)
hold off
end

%% Linear Regression
function [cVal, r2Val, pVal, rmseVal, bVal] = linearRegressionMatlab(time, emgData)

tbl = table(time',emgData','VariableNames',{'TIME','EMGData'});
lm = fitlm(tbl, 'EMGData~TIME','RobustOpts','off'); 
coeff = lm.Coefficients;
cVal   = coeff{'TIME', 'Estimate'};
r2Val  = lm.Rsquared.Ordinary;
pVal   = coeff{'TIME', 'pValue'};
rmseVal= lm.RMSE;
bVal   = coeff{'(Intercept)', 'Estimate'};
end

%% Group Bar plot with error bar 
function plotErrorbarGroup(matMean, matStd)

h = bar(matMean,'hist'); 
hold on
[~, numbars] = size(matMean); 
xdata = get(h,'XData'); 
centerX = cellfun(@(x)(x(1,:)+x(3,:))/2,xdata,'UniformOutput', false);
for i = 1:numbars
    errorbar(centerX{i,:}, matMean(:,i), matStd(:,i), ...
        'linestyle', 'none','LineWidth',2);
end
yline(0,'k');
hold off
xlim([0.5 size(matMean, 1)+0.5]);
set(gca, 'xtick', 1:1:size(matMean, 1));
end