%% import data
clear; clc;

% import benthowave data
dataDir = '20210624 sensor measurement';
experiment = 'Circuit';
numRuns = 2;
circuitData = cell(numRuns,1);

runSkip = 134;
runLabels = 5932:runSkip:6066;

% import all benthowave data into the cell array
for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    circuitData{run} = importfile(dataPath);
end

% import spice
spice = readtable('charge_amp_LT1792_three_stages_calibration.txt');

%%

figure(1);

for run = 1
    data = circuitData{run};
    vout = 10.^(data.Ch2MagdB/20);
    vin = 10.^(data.Ch3MagdB/20);
    
    voutPhaseDeg = 360*data.Ch2Phasecyc;
    vinPhaseDeg = 360*data.Ch3Phasecyc;
    
    subplot(211);
    hold on;
    plot(data.FreqHz, vout./vin, 'r.-');
    plot(spice.Freq, 10.^(spice.Mag/20), 'k-');
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    ylabel('mag [V/V]');

    subplot(212);
    hold on;
    plot(data.FreqHz, voutPhaseDeg-vinPhaseDeg-360, 'r.-');
    plot(spice.Freq, unwrap(spice.Phase), 'k-');
    set(gca,'XScale','log');
    ylabel('phase [deg]');
end

xlabel('frequency [Hz]');
legend('measured','spice');
improvePlot();