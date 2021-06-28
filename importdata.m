%% import data
clear; clc;

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

% data summary:
% proceeds by 2 mm increments
% 1. above water
% 2. at water line
% 3. 2 mm down
% 4. 4 mm down
% ...
% 11. 18 mm down