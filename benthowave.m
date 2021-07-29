importdata();

%% Helmholtz equation
shakeFreq = logspace(1,4);
w = 2*pi*shakeFreq; % angular frequency
h = 10e-3; % depth of sensor
L = 20e-3; % length of vial
g = 9.8;
a = 1;
x0 = a./(w.^2);

c = 2000; % speed of sound in water
rho = 1000; % mass density of water

k = shakeFreq/c;

pressure_acceleration = rho*a.*sin(k*h)./(k.*cos(k*L))+rho*g*x0;

loglog(shakeFreq, abs(pressure_acceleration), '.-');


%% process and plot data
% summary:
% fig. 1: frequency waterfall
% fig. 2: acceleration and predicted pressure at 10 mm depth
% fig. 3: pressure vs depth at 1 kHz

close all;

% frequency waterfall
figure(1);
hold on;

numRuns = 11;
for run = 1:2:numRuns
    data = benthowaveData{run};

    % convert voltage to pressure
    p2v = 10^(-223.3/20)*1e6; % BII-7181 [V/Pa]
    preampGain = 10^(60/20);  % BII-1092

    benthowaveVoltage = 10.^(data.Ch1MagdB/20); % [V}
    pressure = benthowaveVoltage/p2v/preampGain; % [Pa]

    if run == 1
        plot(data.FreqHz, pressure, 'k.-');
    else
        plot(data.FreqHz, pressure, '.-');
    end
end

legend('above water','0 mm','4 mm', '8 mm', '12 mm', '16 mm');

set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('frequency [Hz]');
ylabel('pressure [Pa]');
improvePlot();


% predicted vs measured pressure at 12 mm
figure(2);
hold on;

for run = 1
    data = benthowaveData{run};

    ldvVoltage = 10.^(data.Ch0MagdB/20); % [V]

    vel2volt = 2e-3; % [mm/s/V]

    velocity = vel2volt*ldvVoltage; % [m/s]
    position = 1./(2*pi*data.FreqHz).*velocity; % [m]
    acceleration = 2*pi*data.FreqHz.*velocity; % [m/s^2]

    rho = 1000; % [kg/m^3]
    g = 9.8; % [m/s^2]

    h = 12e-3; % depth [m]

    pressureTheoritical = abs(rho*position.*(g-(2*pi*data.FreqHz).^2*h));
    % Helmholtz equation
    x0 = 1e-3; % displacement amplitude of the shaker
    shakeFreq = data.FreqHz;
    w = 2*pi*shakeFreq; % angular frequency
    L = 26e-3; % length of vial

    c = 2000; % speed of sound in water
    rho = 1000; % mass density of water

    k = shakeFreq/c;

    pressure_acceleration = rho*(g-w.^2.*sin(k*h)./(k.*cos(k*L)));
    pressureHelmholtz = abs(position.*pressure_acceleration);
    
    plot(data.FreqHz, pressureTheoritical, 'r.-');
    plot(data.FreqHz, pressureHelmholtz, 'g.-');
end

% plot for 12 mm depth
for run = 8
    data = benthowaveData{run};

    % convert voltage to pressure
    p2v = 10^(-223.3/20)*1e6; % BII-7181 [V/Pa]
    preampGain = 10^(60/20);  % BII-1092

    benthowaveVoltage = 10.^(data.Ch1MagdB/20); % [V}
    pressure = benthowaveVoltage/p2v/preampGain; % [Pa]

    plot(data.FreqHz, pressure, 'k.-');
end

set(gca,'XScale','log');
set(gca,'YScale','log');

xlabel('frequency [Hz]');
ylabel('pressure [Pa]');
legend('LDV h=12 mm', 'Helmholtz', 'Benthowave h=12 mm');
improvePlot();


% pressure vs depth
figure(3);
hold on;

depths = 0:2:18;
pressures = [];

for run = 2:numRuns
    data = benthowaveData{run};

    idx = find(1000 == data.FreqHz);

    % convert voltage to pressure
    p2v = 10^(-223.3/20)*1e6; % BII-7181 [V/Pa]
    preampGain = 10^(60/20);  % BII-1092

    benthowaveVoltage = 10.^(data.Ch1MagdB/20); % [V}
    pressure = benthowaveVoltage/p2v/preampGain; % [Pa]
    pressures = [pressures pressure(idx)];
end

% measured pressure
plot(depths,pressures, '.-');

% calculated pressure
pressureVsDepth = abs(rho*position(idx).*(g-(2*pi*1000)^2.*(depths*1e-3)));
plot(depths,pressureVsDepth,'r.-');

legend('measured 1 kHz','theoritical 1 kHz');

xlabel('depth [mm]');
ylabel('pressure [Pa]');
improvePlot();