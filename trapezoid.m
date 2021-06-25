%% import data
clear; clc;

dataDir = '20210624 sensor measurement';
experiment = 'Trapezoid';
numRuns = 11;
trapezoidData = cell(numRuns,1);

runSkip = 134;
runLabels = 3905:runSkip:5245;

% import all trapezoid data into the cell array
for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    trapezoidData{run} = importfile(dataPath);
end

% import benthowave data
dataDir = '20210624 sensor measurement';
experiment = 'Benthowave';
numRuns = 11;
benthowaveData = cell(numRuns,1);

runSkip = 134;
runLabels = 2431:runSkip:3771;

% import all benthowave data into the cell array
for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    benthowaveData{run} = importfile(dataPath);
end

%% process and plot data
% summary:
% fig. 1: voltage in vessel
close all;

figure(1);
hold on;

frequency = trapezoidData{1}.FreqHz;

%%%%%%%%%%%%%%%%%%%%%%
% baselines
%%%%%%%%%%%%%%%%%%%%%%
dataAboveWater = trapezoidData{1};
dataAtWaterLine = trapezoidData{2};

% need to convert to volts
plot(frequency, 10.^(dataAboveWater.Ch2MagdB/20), '.-');
plot(frequency, 10.^(dataAtWaterLine.Ch2MagdB/20), '.-');


%%%%%%%%%%%%%%%%%%%%%%
% 5 mm depth
%%%%%%%%%%%%%%%%%%%%%%
N = 3; % number of samples
data5mm = zeros(length(dataAboveWater.Ch2MagdB), N);
for run = 3:5
   data5mm(:,run-2) = 10.^(trapezoidData{run}.Ch2MagdB/20); % [V]
end

% compute mean and standard deviation for each row
data5mmStdDev = std(data5mm,[],2);
data5mmMean = mean(data5mm,2);

% plot 5 mm depth
%plot(frequency, data5mmMean, '.-');
errorbar(frequency,data5mmMean,data5mmStdDev);


%%%%%%%%%%%%%%%%%%%%%%
% 10 mm depth
%%%%%%%%%%%%%%%%%%%%%%
data10mm = zeros(length(dataAboveWater.Ch2MagdB), N);
for run = 6:8
   data10mm(:,run-5) = 10.^(trapezoidData{run}.Ch2MagdB/20);
end

% compute mean and standard deviation for each row
data10mmStdDev = std(data10mm,[],2);
data10mmMean = mean(data10mm,2);

% plot 10 mm depth
%plot(frequency, data10mmMean, '.-');
errorbar(frequency,data10mmMean,data10mmStdDev);


%%%%%%%%%%%%%%%%%%%%%%
% 15 mm depth
%%%%%%%%%%%%%%%%%%%%%%
data15mm = zeros(length(dataAboveWater.Ch2MagdB), N);
for run = 9:11
   data15mm(:,run-8) = 10.^(trapezoidData{run}.Ch2MagdB/20);
end

% compute mean and standard deviation for each row
data15mmStdDev = std(data15mm,[],2);
data15mmMean = mean(data15mm,2);

% plot 10 mm depth
%plot(frequency, data15mmMean, '.-');
errorbar(frequency,data15mmMean,data15mmStdDev);

set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('frequency [Hz]');
ylabel('voltage [V]');

legend('above water','water line','5 mm N=3','10 mm N=3','15 mm N=3');

improvePlot();


figure(2);
hold on;
%%%%%%%%%%%%%%%%%%%%%%
% 10 mm depth
%%%%%%%%%%%%%%%%%%%%%%
data10mm = zeros(length(dataAboveWater.Ch2MagdB), N);
data10mmPhaseDeg = zeros(length(dataAboveWater.Ch2Phasecyc), N);
for run = 6:8
   data10mm(:,run-5) = 10.^(trapezoidData{run}.Ch2MagdB/20);
   data10mmPhaseDeg(:,run-5) = 360*trapezoidData{run}.Ch2Phasecyc;
end

% compute mean and standard deviation for each row
data10mmStdDev = std(data10mm,[],2);
data10mmMean = mean(data10mm,2);
data10mmPhaseDegMean = mean(data10mmPhaseDeg,2);

% plot 10 mm depth
plot(frequency, data10mmMean, '.-');
%errorbar(frequency,data10mmMean,data10mmStdDev);

% plot PVDF minus pickup
plot(frequency, data10mmMean - 10.^(dataAboveWater.Ch2MagdB/20), 'g.-');

ylabel('voltage [V]');
improvePlot();

yyaxis right
% plot for 10 mm depth
for run = 7
    data = benthowaveData{run};

    % convert voltage to pressure
    p2v = 10^(-223.3/20)*1e6; % BII-7181 [V/Pa]
    preampGain = 10^(60/20);  % BII-1092

    benthowaveVoltage = 10.^(data.Ch1MagdB/20); % [V}
    pressure = benthowaveVoltage/p2v/preampGain; % [Pa]
    pressurePhaseDeg = 360*data.Ch1Phasecyc;

    plot(data.FreqHz, pressure, '.-');
end

set(gca,'XScale','log');
set(gca,'YScale','log');

ylabel('pressure [Pa]');
xlabel('frequency [Hz]');

improvePlot();


% pressure to voltage transfer function
figure(3);
hold on;
subplot(211);
plot(frequency, data10mmMean./pressure, '.-');

set(gca,'XScale','log');
set(gca,'YScale','log');

xlim([100 4400]);
ylabel('mag [V/Pa]');
xlabel('frequency [Hz]');

subplot(212);
plot(frequency, data10mmPhaseDegMean - pressurePhaseDeg, '.-');

xlim([100 4400]);

set(gca,'XScale','log');
ylabel('phase [deg]');
xlabel('frequency [Hz]');

improvePlot();

