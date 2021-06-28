%% import data
importdata();

%% analytical circuit

f = wormData{1}.FreqHz;

s = tf('s');

R0 = 10e6;

R1 = 34e6;
C1 = 10e-12;

f1 = 1/(2*pi*R1*C1);

R2 = 10e3;
C2 = 100e-9;

f2 = 1/(2*pi*R2*C2);

R3 = 100e3;
C3 = 50e-12;

f3 = 1/(2*pi*R3*C3);

R4 = 1e3;
R5 = 100e3;

% first stage
Z1 = R1/(s*C1)/(R1+1/(s*C1));
S1 = -s*Z1; % transimpedance

% second stage
Z2 = R2 + 1/(s*C2);
Z3 = R3/(s*C3)/(R3+1/(s*C3));

S2 = -Z3/Z2;

% final stage
S3 = -R5/R4;

[m, p] = bode(S1*S2*S3, 2*pi*f);


%% worm analytical
COMSOL_freq_top_electrode_import();

w = 0.5e-3;
l = 15e-3;
d31_prime = 22e-12;
d32_prime = 3e-12;
d33_prime = -30e-12;

E = 5.7e9;
% add damping by making E complex
zeta = [0 0.1 0.5 1];
figure;
subplot(2,1,1);
ylabel('Magnitude [C/Pa]')
for index = 1:length(zeta)
    E = E + zeta(index)*E*1i;
    nu = 0.35;
    rho = 1780;

    c = sqrt(E/rho);

    P0 = 1; % Pascals

    piezo_angle = 0;
    d31 = cos(piezo_angle).^2*d31_prime+sin(piezo_angle).^2*d32_prime;
    d32 = sin(piezo_angle).^2*d31_prime+cos(piezo_angle).^2*d32_prime;
    d33 = d33_prime*ones(size(piezo_angle));

    freq = logspace(log10(200),6,600);
    ang_freq = 2*pi*freq_comsol';
    k = ang_freq/c;

    F11 = -P0*w*((1-2*nu)*tan(k*l)./k+2*nu*l);
    F22 = -P0*w*l*ones(size(ang_freq));
    F33 = F22;

    Q31 = d31*F11;
    Q32 = d32*F22;
    Q33 = d33*F33;

    Q = Q31+Q32+Q33;

    hearing_range = freq_comsol < 20e3;

    %cost_func = 'NRMSE';
    %fit = goodnessOfFit(abs((Q(hearing_range))'),abs(Q_comsol((hearing_range))),cost_func);
    %error = (1-fit)*100;

    subplot(2,1,1);
    hold on;
    plot(freq_comsol,abs(Q),'-');
    xlim([200 max(freq_comsol)]);
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    
    subplot(2,1,2);
    hold on;
    plot(freq_comsol,180/pi*angle(Q),'-');
    set(gca, 'YScale', 'linear')
    set(gca, 'XScale', 'log')
    xlabel('Frequency [Hz]')
    ylabel('Phase [^\circ]')
end
hold off;
subplot(2,1,1);
plot(freq_comsol,abs(Q_comsol),'k.');
legend([ '\zeta = ' + string(zeta) 'COMSOL']);
    
improvePlot();

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
magicNum = 1e-14;
plot(frequency, magicNum*m(:), 'k.-');

% plot new low freq data
dataLowFreq = readtable('low_freq_PVDF.csv');
plot(dataLowFreq.f, dataLowFreq.m./(pressure(1:length(dataLowFreq.m))), 'b.-');

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
