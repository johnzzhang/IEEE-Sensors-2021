%% import data
clear; clc;

dataDir = '20210624 sensor measurement';
experiment = 'Worm';
numRuns = 15;
wormData = cell(numRuns,1);

runSkip = 134;
runLabels = 421:runSkip:2297;

% import all trapezoid data into the cell array
for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    wormData{run} = importfile(dataPath);
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

% import spice
spice = readtable('charge_amp_LT1792_three_stages.txt');

%% process and plot

% summary:
% fig. 1: voltage in vessel
close all;

figure(1);
hold on;

frequency = wormData{1}.FreqHz;

%%%%%%%%%%%%%%%%%%%%%%
% baselines
%%%%%%%%%%%%%%%%%%%%%%
dataAtWaterLine = wormData{1};

%plot(frequency, 10.^(dataAtWaterLine.Ch2MagdB/20), '.-');

%%%%%%%%%%%%%%%%%%%%%%
% 10 mm depth
%%%%%%%%%%%%%%%%%%%%%%
N = 3;
data10mm = zeros(length(dataAtWaterLine.Ch2MagdB), N);
data10mmPhaseDeg = zeros(length(dataAtWaterLine.Ch2Phasecyc), N);
for run = 6:8
   data10mm(:,run-5) = 10.^(wormData{run}.Ch2MagdB/20);
   data10mmPhaseDeg(:,run-5) = 360*wormData{run}.Ch2Phasecyc;
end

% compute mean and standard deviation for each row
data10mmStdDev = std(data10mm,[],2);
data10mmMean = mean(data10mm,2);
data10mmPhaseDegMean = mean(data10mmPhaseDeg,2);

% plot 10 mm depth
%plot(frequency, data10mmMean, '.-');
errorbar(frequency,data10mmMean,data10mmStdDev);

set(gca,'XScale','log');
set(gca,'YScale','log');

xlim([100 10000]);

yyaxis right

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
xlim([100 10000]);

improvePlot();

% pressure to voltage transfer function
figure(2);
subplot(211);
hold on;
plot(frequency, data10mmMean./pressure, 'r.-');

% plot charge amp and 
magicNum = 2e-15;
plot(spice.Freq, magicNum*10.^(spice.Mag/20), 'k.-');

set(gca,'XScale','log');
set(gca,'YScale','log');

ylabel('mag [V/Pa]');
xlabel('frequency [Hz]');

subplot(212);
plot(frequency, data10mmPhaseDegMean - pressurePhaseDeg, 'r.-');

set(gca,'XScale','log');
ylabel('phase [deg]');
xlabel('frequency [Hz]');

improvePlot();
