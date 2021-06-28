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

%% analytical circuit

f = logspace(1,log10(100e3),20);

s = tf('s');

R0 = 10e6;

R1 = 34e6;
C1 = 10e-12;

f1 = 1/(2*pi*R1*C1);

R2 = 10e3;
C2 = 100e-9;

f2 = 1/(2*pi*R2*C2);

R3 = 100e3;
%C3 = 100e-12;
C3 = 50e-12;

f3 = 1/(2*pi*R3*C3);

R4 = 1e3;
R5 = 100e3;

% first stage
Z1 = R1/(s*C1)/(R1+1/(s*C1));
S1 = -Z1/R0;

% second stage
Z2 = R2 + 1/(s*C2);
Z3 = R3/(s*C3)/(R3+1/(s*C3));

S2 = -Z3/Z2;

% final stage
S3 = -R5/R4;

[m, p] = bode(S1*S2*S3, 2*pi*f);

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
    plot(data.FreqHz, vout./vin, 'rx');
    plot(f, m(:), 'k-');
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    ylabel('mag [V/V]');

    subplot(212);
    hold on;
    plot(data.FreqHz, voutPhaseDeg-vinPhaseDeg, 'rx');
    plot(f, p(:), 'k-');
    set(gca,'XScale','log');
    ylabel('phase [deg]');
end

xlabel('frequency [Hz]');
legend('measured','spice');
improvePlot();