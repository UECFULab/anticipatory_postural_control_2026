function analysisMotion()
clear all
close all

matDir = 'data/';
figDir = 'figs/';
if ~exist(figDir, 'dir')
    mkdir(figDir);
end
saveIndFigs = 0;

extractDur = 5; % Duration (s) for analysis with Floor tilt start as 0 (s)
sList = 1:10;
% sList = 1:2;

diffNonSList = cell(1,length(sList));
diffConSList = cell(1,length(sList));
diffNonSFloorList = cell(1,length(sList));
diffConSFloorList = cell(1,length(sList));
%
diffAnkleAngNonSList = cell(1,length(sList));
diffAnkleAngConSList = cell(1,length(sList));
diffHipAngNonSList = cell(1,length(sList));
diffHipAngConSList = cell(1,length(sList));
%
COMMeanSubNon = zeros(length(sList), 1);
COMMeanSubCon = zeros(length(sList), 1);
angAnkleMeanSubNon = zeros(length(sList), 1);
angAnkleMeanSubCon = zeros(length(sList), 1);
angHipMeanSubNon = zeros(length(sList), 1);
angHipMeanSubCon = zeros(length(sList), 1);
%
COMMeanSub = zeros(length(sList), 1);
angAnkleMeanSub = zeros(length(sList), 1);
angHipMeanSub = zeros(length(sList), 1);
%
COMMeanConAllSub = [];
COMMeanAllSub = [];
angAnkleConAllSub = [];
angHipConAllSub   = [];
%
for sIndx = sList
    Name = ['subject' num2str(sIndx)]; % Subject Number
    if(strcmp(Name, 'subject1'))
        SkipData = 2;
        HipAngRange = [0.2,0.4];
    elseif(strcmp(Name, 'subject2'))
        SkipData = [];
        HipAngRange = [0.1,0.3];
    elseif(strcmp(Name, 'subject3'))
        SkipData = [1, 57];
        HipAngRange = [0.1,0.3];
    elseif(strcmp(Name, 'subject4'))
        SkipData = [];
        HipAngRange = [0.1,0.3];
    elseif(strcmp(Name, 'subject5'))
        SkipData = 10;
        HipAngRange = [0.1,0.3];
    elseif(strcmp(Name, 'subject6'))
        SkipData = [16 20];
        HipAngRange = [0.0,0.2];
    elseif(strcmp(Name, 'subject7'))
        SkipData = 1;
        HipAngRange = [0.0,0.2];
    elseif(strcmp(Name, 'subject8'))
        SkipData = [];
        HipAngRange = [0.0,0.2];
    elseif(strcmp(Name, 'subject9'))
        SkipData = 44;
        HipAngRange = [0.0,0.2];
    elseif(strcmp(Name, 'subject10'))
        SkipData = 1;
        HipAngRange = [0.0,0.2];
    else
        SkipData = [];
    end
    %-- Read measured motion data --% 
    matFileName = [matDir, 'MotionData_',Name,'.xlsx'];
    comX        = readmatrix(matFileName,  'sheet', 'comX');
    angAnkle    = readmatrix(matFileName,  'sheet', 'angAnkle');
    angHip      = readmatrix(matFileName,  'sheet', 'angHip');
    angForcePlate = readmatrix(matFileName,'sheet', 'angForcePlate');
    posAnkleX   = readmatrix(matFileName,  'sheet', 'posAnkleX');
    tiltEndIndx = readmatrix(matFileName,  'sheet', 'tiltEndIndx'); % measured freqeuency: 300Hz
    %
    Soundstart  = readmatrix(matFileName,  'sheet', 'Soundstart');  % measured freqeuency: 1000Hz
    Soundend    = readmatrix(matFileName,  'sheet', 'Soundend');    % measured freqeuency: 1000Hz
    %--%

    nonIndx = 11:30;
    conIndx = 41:60;
    nonIndx = nonIndx(~ismember(nonIndx, SkipData)); % Remove failed trial (SkipData)
    conIndx = conIndx(~ismember(conIndx, SkipData)); % Remove failed trial (SkipData)

    %- Set ankle position as 0 mm (for COM position) -%
    comX  = comX-posAnkleX;
    %--%

    %
    COMTime = linspace(-extractDur, extractDur, 2*extractDur*300);
    Soundstartmean = mean(Soundstart(conIndx),2)/1000;
    Soundendmean   = mean(Soundend(conIndx),2)/1000;
    SoundstartSD   = std(Soundstart(conIndx),[], 2)/1000;
    SoundendSD     = std(Soundend(conIndx),[], 2)/1000;
    tiltEndMean    = mean(tiltEndIndx(conIndx),  2)/300;

    disp(['Cue start time (mean(SD)) : ' num2str(Soundstartmean, '%.2f') '(' num2str(SoundstartSD, '%.2f') ') (s), Cue end time (mean) : ' num2str(Soundendmean) ' (s), Tilt end time (mean) : ' num2str(tiltEndMean) ' (s)']);

    % 
    figure;
    title(['Ankle_' Name])
    plotResultInd(COMTime, angAnkle(:,nonIndx), angAnkle(:,conIndx), Soundstartmean, Soundendmean, SoundstartSD, SoundendSD, tiltEndMean, [-0.2, -0.0], -0.2:0.1:0.0, false, extractDur)
    if(saveIndFigs)
        saveas(gcf, [figDir 'Ankle_' Name],'emf')
    end
    % Plot angHipRE
    figure;
    title(['Hip_' Name])
    plotResultInd(COMTime, angHip(:,nonIndx), angHip(:,conIndx), Soundstartmean, Soundendmean, SoundstartSD, SoundendSD, tiltEndMean, HipAngRange, 0.0:0.1:0.4, false, extractDur)
    if(saveIndFigs)
        saveas(gcf, [figDir 'Hip_' Name],'emf')
    end
    % Plot angHipRE
    figure;
    title(['COM_' Name])
    plotResultInd(COMTime, angForcePlate(:,nonIndx), angForcePlate(:,conIndx), Soundstartmean, Soundendmean, SoundstartSD, SoundendSD, tiltEndMean, [], [], true, extractDur)
    title(['ForcePlate Subject: ' num2str(sIndx)]);
    %
    figure;
    plotResultInd(COMTime, comX(:,nonIndx), comX(:,conIndx), Soundstartmean, Soundendmean, SoundstartSD, SoundendSD, tiltEndMean, [0 100], 0:50:100, false, extractDur)    % title(['COM-Ankle X Subject: ' num2str(sIndx)]);
    if(saveIndFigs)
        saveas(gcf, [figDir 'COM_' Name],'emf')
    end
    comXRENonS = comX(:,nonIndx);
    comXREConS = comX(:,conIndx);
    
    [~,soundStartIndx] = min(abs(COMTime - Soundstartmean));
    [~,floorStartIndx] = min(abs(COMTime - 0));
    [~,floorStopIndx]  = min(abs(1-COMTime));

    diffNonSList{sIndx} = comXRENonS(floorStartIndx, :) -comXRENonS(soundStartIndx, :);
    diffConSList{sIndx} = comXREConS(floorStartIndx, :) -comXREConS(soundStartIndx, :);
    diffNonSFloorList{sIndx} = comXRENonS(floorStopIndx, :) -comXRENonS(soundStartIndx, :);
    diffConSFloorList{sIndx} = comXREConS(floorStopIndx, :) -comXREConS(soundStartIndx, :);
    %
    diffAnkleAngNonSList{sIndx} = angAnkle(floorStartIndx, nonIndx) -angAnkle(soundStartIndx, nonIndx);
    diffAnkleAngConSList{sIndx} = angAnkle(floorStartIndx, conIndx) -angAnkle(soundStartIndx, conIndx);
    diffHipAngNonSList{sIndx} = angHip(floorStartIndx, nonIndx) -angHip(soundStartIndx, nonIndx);
    diffHipAngConSList{sIndx} = angHip(floorStartIndx, conIndx) -angHip(soundStartIndx, conIndx);
    %
    angAnkleMeanSubNon(sIndx) = mean(mean(angAnkle(1:soundStartIndx,nonIndx), 1));
    angAnkleMeanSubCon(sIndx) = mean(mean(angAnkle(1:soundStartIndx,conIndx), 1));
    angHipMeanSubNon(sIndx) = mean(mean(angHip(1:soundStartIndx,nonIndx), 1));
    angHipMeanSubCon(sIndx) = mean(mean(angHip(1:soundStartIndx,conIndx), 1));
    COMMeanSubNon(sIndx) = mean(mean(comX(1:soundStartIndx,nonIndx), 1));
    COMMeanSubCon(sIndx) = mean(mean(comX(1:soundStartIndx,conIndx), 1));
    %
    angAnkleMeanSub(sIndx) = mean(mean(angAnkle(1:soundStartIndx,:), 1));
    angHipMeanSub(sIndx)   = mean(mean(angHip(1:soundStartIndx,:), 1));
    COMMeanSub(sIndx)      = mean(mean(comX(1:soundStartIndx,:), 1));
    %
    COMMeanConAllSub  = [COMMeanConAllSub, mean(comX(1:soundStartIndx,conIndx), 1)];
    COMMeanAllSub     = [COMMeanAllSub, mean(comX(1:soundStartIndx,:), 1)];%
    angAnkleConAllSub = [angAnkleConAllSub, mean(angAnkle(1:soundStartIndx,conIndx), 1)];
    angHipConAllSub   = [angHipConAllSub, mean(angHip(1:soundStartIndx,conIndx), 1)];
end

%-- Bar plot: difference between CS-FS --%
figure;
subplot(3,1,1)
plotErrorbarGroup([cellfun(@(x) mean(x,2), diffNonSList)', cellfun(@(x) mean(x,2), diffConSList)'], [cellfun(@(x) std(x,0,2), diffNonSList)', cellfun(@(x) std(x,0,2), diffConSList)'], [])
ylim([-10,25])
set(gca, 'ytick', -10:10:30);
ylabel('COM');
%
subplot(3,1,2)
plotErrorbarGroup([cellfun(@(x) mean(x,2), diffHipAngNonSList)', cellfun(@(x) mean(x,2), diffHipAngConSList)'], [cellfun(@(x) std(x,0,2), diffHipAngNonSList)', cellfun(@(x) std(x,0,2), diffHipAngConSList)'], [])
ylim([-0.02,0.02])
set(gca, 'ytick', -0.02:0.01:0.02);
ylabel('Hip');
%
subplot(3,1,3)
plotErrorbarGroup([cellfun(@(x) mean(x,2), diffAnkleAngNonSList)', cellfun(@(x) mean(x,2), diffAnkleAngConSList)'], [cellfun(@(x) std(x,0,2), diffAnkleAngNonSList)', cellfun(@(x) std(x,0,2), diffAnkleAngConSList)'], [])
ylim([-0.02,0.02])
set(gca, 'ytick', -0.02:0.01:0.02);
ylabel('Ankle');
%
saveas(gcf, [figDir 'COMAnkleHipdiff'],'emf')
%--% 

%-- t-Test: difference between CS-FS --%
%- COM -%
tResults = zeros(length(sList),3);
for sIndx = sList
    [~,tResults(sIndx,1),~,states] = ttest2(diffNonSList{sIndx}, diffConSList{sIndx});
    tResults(sIndx,2)  = states.tstat;
    tResults(sIndx,3) = states.df;
end
tResultCOMTable = array2table(tResults,'VariableNames',{'p-value','t-value','df'});
writetable(tResultCOMTable, 'tResultCOM.csv');
%- Hip -%
tResults = zeros(length(sList),3);
for sIndx = sList
    [~,tResults(sIndx,1),~,states] = ttest2(diffHipAngNonSList{sIndx}, diffHipAngConSList{sIndx});
    tResults(sIndx,2)  = states.tstat;
    tResults(sIndx,3) = states.df;
end
tResultHipTable = array2table(tResults,'VariableNames',{'p-value','t-value','df'});
writetable(tResultHipTable, 'tResultHip.csv');
%- Ankle -%
tResults = zeros(length(sList),3);
for sIndx = sList
    [~,tResults(sIndx,1),~,states] = ttest2(diffAnkleAngNonSList{sIndx}, diffAnkleAngConSList{sIndx});
    tResults(sIndx,2)  = states.tstat;
    tResults(sIndx,3) = states.df;
end
tResultAnkTable = array2table(tResults,'VariableNames',{'p-value','t-value','df'});
writetable(tResultAnkTable, 'tResultAnk.csv');
%--%

disp(['COM no cue:', num2str(mean(COMMeanSubNon)), ', with cue: ', num2str(mean(COMMeanSubCon))]);
disp(['Ankle no cue:', num2str(mean(angAnkleMeanSubNon)), ', with cue: ', num2str(mean(angAnkleMeanSubCon))]);
disp(['Hip no cue:', num2str(mean(angHipMeanSubNon)), ', with cue: ', num2str(mean(angHipMeanSubCon))]);
%
disp(['COM Mean:', num2str(mean(COMMeanSub))]);
disp(['Ankle Mean:', num2str(mean(angAnkleMeanSub))]);
disp(['Hip Mean:', num2str(mean(angHipMeanSub))]);
%
disp(['COM at CS for all subjects (Mean(SD)) : ' num2str(mean(COMMeanConAllSub), '%.1f') '(' num2str(std(COMMeanConAllSub), '%.1f') ') (mm) ']);
disp(['Ankle angle at CS for all subjects (Mean(SD)) : ' num2str(mean(angAnkleConAllSub), '%.2f') '(' num2str(std(angAnkleConAllSub), '%.2f') ') (rad) ']);
disp(['Hip angle at CS for all subjects (Mean(SD)) : ' num2str(mean(angHipConAllSub), '%.2f') '(' num2str(std(angHipConAllSub), '%.2f') ') (rad) ']);
end

%% Group Bar plot with error bar 
function plotErrorbarGroup(matMean, matStd, xName)

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
set(gca, 'XTickLabel', xName);
xlim([0.5 size(matMean, 1)+0.5]);
set(gca, 'xtick', 1:1:size(matMean, 1));
end

%% Plot individual data
function plotResultInd(time, matN, matC, Soundstartmean, Soundendmean,  SoundstartSD, SoundendSD, tiltEndMean, yRange, ytickVal ,showAll, extractDur)

matNM  = mean(matN,2);
matNSD = std(matN,0,2);
matCM  = mean(matC,2);
matCSD = std(matC,0,2);

subplot(2,1,1)
if(showAll)
    plot(time, matN, '-', 'LineWidth',1)
else
    plot(time, matNM, 'k-', 'LineWidth',1)
    hold on;
    plot(time, matNM+matNSD, 'c-',  'LineWidth',1)
    plot(time, matNM-matNSD, 'c-',  'LineWidth',1)
end
xline(0,'r');
xline(tiltEndMean,'r-');
hold off
xlim([-extractDur, extractDur])
if(~isempty(yRange))
    ylim(yRange)
    yticks(ytickVal)
end
%
subplot(2,1,2)
if(showAll)
    plot(time, matC, '-', 'LineWidth',1)
else
    plot(time, matCM, 'k-', 'LineWidth',1)
    hold on;
    plot(time, matCM+matCSD, 'c-',  'LineWidth',1)
    plot(time, matCM-matCSD, 'c-',  'LineWidth',1)
end
xline(0,'r');
xline(tiltEndMean,'r-');
xline(Soundstartmean,'b-');
xline(Soundstartmean+SoundstartSD,'c-');
xline(Soundstartmean-SoundstartSD,'c-');
xline(Soundendmean,'b-');
xline(Soundendmean+SoundendSD,'c-');
xline(Soundendmean-SoundendSD,'c-');
hold off
xlim([-extractDur, extractDur])
if(~isempty(yRange))
    ylim(yRange)
    yticks(ytickVal)
end
end
