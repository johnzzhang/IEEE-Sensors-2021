close all; clc; clear;

importdata();

%% process worm data
airDecoyNoise = shieldedWormData{1};
waterDecoyUnshielded = shieldedWormData{2};
waterDecoyShielded = shieldedWormData{3};
waterDecoyNoise = shieldedWormData{4};
waterRealNoise = shieldedWormData{5};
waterRealShielded = shieldedWormData{6};

%% compare noise inside and outside water
figure(1);

hold on;

plot(airDecoyNoise.FreqHz, airDecoyNoise.Ch2MagdB, 'k.-');
plot(waterDecoyNoise.FreqHz, waterDecoyNoise.Ch2MagdB, 'r.-');
plot(waterRealNoise.FreqHz, waterRealNoise.Ch2MagdB, 'b.-');

set(gca,'XScale','log');

legend('air decoy','water decoy','water real');

ylabel('charge amp output [dB V]');
xlabel('frequency [Hz]');

xlim([10, 20e3]);

improvePlot();

%% compare coupling shoot through

figure(2);

hold on;

plot(waterDecoyUnshielded.FreqHz, waterDecoyUnshielded.Ch2MagdB, 'k.-');
plot(waterDecoyShielded.FreqHz, waterDecoyShielded.Ch2MagdB, 'r.-');
plot(waterRealShielded.FreqHz, waterRealShielded.Ch2MagdB, 'b.-');

set(gca,'XScale','log');

legend('decoy unshielded','decoy shielded','real shielded');

ylabel('charge amp output [dB V]');
xlabel('frequency [Hz]');

xlim([200, 12.8e3]);

improvePlot();

%% compare acceleration and charge amp output
figure(3);

hold on;

plot(waterRealShielded.FreqHz, waterRealShielded.Ch2MagdB, 'b.-');
plot(waterRealShielded.FreqHz, waterRealShielded.Ch1MagdB, 'r.-');

set(gca,'XScale','log');

ylabel('[dB V]');
xlabel('frequency [Hz]');

legend('charge amp output', 'accelerometer');

xlim([200, 12.8e3]);

improvePlot();

%% Analytical Model
f = waterRealShielded.FreqHz;

% convert dB to V to Pa
g = 9.8;
volts2accel = g/(10e-3); % [m/s^2]
accelOutput = volts2accel*10.^(waterRealShielded.Ch1MagdB/20);

% create pressure distribution from measured data
pressureDistribution = cell(length(f),1);

rho_water = 1000;
% create interpolated pressure
for i = 1:length(f)
   % gravity + acceleration
   % contribution from gravity is negligible
   pressureDistribution{i} = @(x) rho_water*g*accelOutput(i)/(2*pi*f(i)) + rho_water*x*accelOutput(i);
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

[m, p] = bode(S1*S2*S3, 2*pi*f);


%% plot analytical on measured
figure(4);

hold on;

% analytical
fill([f;flipud(f)],[abs(Qupper);flipud(abs(Qlower))],'k','linestyle','none');
alpha(0.25);

% measured
signal = 10.^(waterRealShielded.Ch2MagdB/20);
plot(f, signal./m(:), 'bx');

set(gca,'XScale','log');
set(gca,'YScale','log');

ylabel('charge output [C]');

set(gca,'XScale','log');
xlabel('frequency [Hz]');
improvePlot();

