close all; clc; clear;

importdata();

%% Analytical Model
f = benthowaveData{1}.FreqHz;
depths = 0:2:18;
depths = 1e-3*depths;

% create pressure distribution from measured data
pressureDistribution = cell(length(f),1);

pressureMatrix = zeros(length(depths), length(f));
pressureScale = 1; % assume linear relation to voltage input to shaker
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
   pressureDistribution{i} = @(x) pressureScale.*interp1(depths,pressureMatrix(:,i),x);
end

% analytical worm
w = 1e-3;
l = 18e-3;
d31 = 3e-12;
d32 = d31;
d33 = -40e-12;

E = 5.6e9;
nu = 0.37;
rho = 1780;

s11 = 1/E;
s12 = -nu/E;
s13 = -nu/E;

c = sqrt(E/rho);

k = (2*pi*f)/c;

h = 18e-3;

% part without pressure varying
Q = zeros(size(f));
for i = 1:length(Q)
    Q = -d31/s11*w*(s11+s12+s13).*pressureDistribution{i}(h).*tan(k*l)./k;
end
Qpressure = zeros(size(Q));
Qupper = zeros(size(Q));

for i = 1:length(Q)
    Qpressure(i) = w*((s12+s13)/s11*d31-d32-d33)*integral(pressureDistribution{i}, 0, h, 'ArrayValued', true);
    Qupper(i) = Q(i) + Qpressure(i);
end

% analytical worm
w = 1e-3;
l = 18e-3;
d31 = 5e-12;
d32 = d31;
d33 = -20e-12;

E = 5.6e9;
nu = 0.37;
rho = 1780;

s11 = 1/E;
s12 = -nu/E;
s13 = -nu/E;

c = sqrt(E/rho);

k = (2*pi*f)/c;

h = 18e-3;

% part without pressure varying
Q = zeros(size(f));
for i = 1:length(Q)
    Q = -d31/s11*w*(s11+s12+s13).*pressureDistribution{i}(h).*tan(k*l)./k;
end
Qpressure = zeros(size(Q));
Qlower = zeros(size(Q));

for i = 1:length(Q)
    Qpressure(i) = w*((s12+s13)/s11*d31-d32-d33)*integral(pressureDistribution{i}, 0, h, 'ArrayValued', true);
    Qlower(i) = Q(i) + Qpressure(i);
end

%% process new worm data
totalData = newWormData{1}; % PVDF in water with pressure
couplingData = newWormData{2}; % PVDF in water no pressure

frequency = totalData.FreqHz;

total = 10.^(totalData.Ch1MagdB/20);
coupling = 10.^(couplingData.Ch1MagdB/20);
signal = total-coupling;

%% analytical circuit
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

[m, p] = bode(S1*S2*S3, 2*pi*frequency);


%% compare shootthrough noise to signal
figure(1);

xlims = [1e3 10e3];

%subplot(211);
hold on;

plot(frequency, total, 'k.-');
plot(frequency, coupling, 'r.-');

set(gca,'XScale','log');
%set(gca,'YScale','log');

legend('total','coupling');

ylabel('charge amp output [V]');
xlabel('frequency [Hz]');

xlim(xlims);
improvePlot();

%% plot analytical on measured
figure(2);

hold on;

% analytical
fill([f;flipud(f)],[abs(Qupper);flipud(abs(Qlower))],'k','linestyle','none');
alpha(0.25);

% measured
plot(frequency, signal./m(:), 'bx');

set(gca,'XScale','log');
set(gca,'YScale','log');

ylabel('charge output [C]');

xlim(xlims);
set(gca,'XScale','log');
xlabel('frequency [Hz]');
improvePlot();

