importdata();

%% Analytical Model
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
s13 = -nu/E;

c = sqrt(E/rho);

k = (2*pi*f)/c;

h = 10e-3;

% part without pressure varying
Q = zeros(size(f));
for i = 1:length(Q)
    Q = -d31/s11*w*(s11+s12+s13).*pressureDistribution{i}(h).*tan(k*l)./k;
end
Qpressure = zeros(size(Q));
Qtotal = zeros(size(Q));

for i = 1:length(Q)
    Qpressure(i) = w*((s12+s13)/s11*d31-d32-d33)*integral(pressureDistribution{i}, 0, h, 'ArrayValued', true);
    Qtotal(i) = Q(i) + Qpressure(i);
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

frequency = wormData{1}.FreqHz;

%%%%%%%%%%%%%%%%%%%%%%
% baselines
%%%%%%%%%%%%%%%%%%%%%%
dataAtWaterLine = wormData{1};

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


for run = 7
    data = benthowaveData{run};

    % convert voltage to pressure
    p2v = 10^(-223.3/20)*1e6; % BII-7181 [V/Pa]
    preampGain = 10^(60/20);  % BII-1092

    benthowaveVoltage = 10.^(data.Ch1MagdB/20); % [V}
    pressure = benthowaveVoltage/p2v/preampGain; % [Pa]
    pressurePhaseDeg = 360*data.Ch1Phasecyc;
end

%% pressure to voltage transfer function
figure(2);

xlims = [200 20e3];

%subplot(211);
hold on;
plot(f, data10mmMean./pressure, 'r.-');
%errorbar(frequency,data10mmMean./pressure,data10mmStdDev./pressure);

% plot voltage from analytical
%plot(f, m(:).*abs(Qpressure), 'm.-');
%plot(f, m(:).*abs(Q), 'g.-');
plot(f, m(:).*abs(Qtotal), 'k.-');

legend('measured','theory');

xlim(xlims);

set(gca,'XScale','log');
set(gca,'YScale','log');

ylabel('charge amp output [V]');
xlabel('frequency [Hz]');

% plotting phase doesn't make sense for frequency response
% subplot(212);
% hold on;
% wormPhase = data10mmPhaseDegMean - pressurePhaseDeg;
% wormPhase(2) = wormPhase(1) + 360; % manual unwrap
% plot(f, unwrap(wormPhase)+360, 'r.-');
% plot(f, p(:), 'k.-');

xlim(xlims);
set(gca,'XScale','log');
%ylabel('phase [deg]');
xlabel('frequency [Hz]');

improvePlot();
