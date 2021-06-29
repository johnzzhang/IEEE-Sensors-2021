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

close all;

%% Analytical
f = benthowaveData{1}.FreqHz;
depths = 0:2:18;
depths = 1e-3*depths;

% create pressure distribution from measured data
pressureDistribution = cell(length(f),1);

pressureMatrix = zeros(length(depths), length(f));
% 0 to 18 [mm]
for run = 1:length(depths)
    data = benthowaveData{run};

    % convert voltage to pressure
    p2v = 10^(-223.3/20)*1e6; % BII-7181 [V/Pa]
    preampGain = 10^(60/20);  % BII-1092

    benthowaveVoltage = 10.^(data.Ch1MagdB/20); % [V}
    pressure = benthowaveVoltage/p2v/preampGain; % [Pa]

    pressureMatrix(run,:) = pressure;
end

% create interpolated pressure
for i = 1:length(f)
   pressureDistribution{i} = @(x) interp1(depths,pressureMatrix(:,i),x);
end

% analytical worm
w = 0.5e-3;
l = 15e-3;
d31 = 3e-12;
d32 = 3e-12;
d33 = -30e-12;

E = 5.6e9;
nu = 0.37;
rho = 1780;

s11 = 1/E;
s12 = -nu/E;
s13 = s12;

c = sqrt(E/rho);

k = (2*pi*f)/c;

h = 10e-3;

% part without pressure varying
Q = -w*d31/s11*(s11+s12+s13).*pressureDistribution{i}(h).*tan(k*l)./k;
Qpressure = zeros(size(Q));

for i = 1:length(Q)
    Qpressure(i) = Q(i) + w*((s12+s13)/s11*d31-d32-d33)*integral(pressureDistribution{i}, 0, 10e-3, 'ArrayValued', true);
end

% analytical circuit
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

%% process and plot

% summary:
% fig. 1: voltage in vessel
% close all;
% 
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
close all;

%% pressure to voltage transfer function
figure(2);

xlims = [200 20e3];

subplot(211);
hold on;
plot(f, data10mmMean./pressure, 'r.-');

% plot voltage from analytical
plot(f, m(:).*abs(Qpressure), 'k.-');

xlim(xlims);

set(gca,'XScale','log');
set(gca,'YScale','log');

ylabel('mag [V/Pa]');
xlabel('frequency [Hz]');

subplot(212);
wormPhase = data10mmPhaseDegMean - pressurePhaseDeg;
wormPhase(2) = wormPhase(1) + 360; % manual unwrap
plot(frequency, unwrap(wormPhase), 'r.-');

xlim(xlims);
set(gca,'XScale','log');
ylabel('phase [deg]');
xlabel('frequency [Hz]');

improvePlot();
