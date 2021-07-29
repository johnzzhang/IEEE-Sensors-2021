%% import data
close all; clear; clc;

%% import all benthowave data into the cell array
dataDir = '20210624 sensor measurement';
experiment = 'Benthowave';
numRuns = 11;
benthowaveData = cell(numRuns,1);

runSkip = 134;
runLabels = 2431:runSkip:3771;

for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    benthowaveData{run} = importfile(dataPath);
end

% data summary:
% proceeds by 2 mm increments
% 1. above water
% 2. at water line
% 3. 2 mm down
% 4. 4 mm down
% ...
% 11. 18 mm down

%% import all worm data into the cell array
dataDir = '20210624 sensor measurement';
experiment = 'Worm';
numRuns = 15;
wormData = cell(numRuns,1);

runSkip = 134;
runLabels = 421:runSkip:2297;


for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    wormData{run} = importfile(dataPath);
end


%% import all noise data into the cell array

dataDir = '20210624 sensor measurement';
experiment = 'Noise';
numRuns = 4;
noiseData = cell(numRuns,1);

runSkip = 134;
runLabels = 5379:runSkip:5781;

for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/' experiment '/f' runNumber '.txt'];
    noiseData{run} = importfile(dataPath);
end

%% import new data for worm
dataDir = '20210715';
numRuns = 2;

runSkip = 37;
runLabels = [253 389];

newWormData = cell(numRuns,1);

for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/f' runNumber '.txt'];
    newWormData{run} = importfile(dataPath);
end

%% import new shielded data
dataDir = '20210727';

runSkip = 37;
runLabels = [10 144 300 301 305 392];

numRuns = numel(runLabels);

shieldedWormData = cell(numRuns,1);

for run = 1:numRuns
    runNumber = num2str(runLabels(run));
    dataPath = ['data/' dataDir '/f' runNumber '.txt'];
    shieldedWormData{run} = importfile(dataPath);
end
