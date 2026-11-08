% Problem 6.2 and 6.3: static and dynamic Kalman filters.
% 6.2: 4-state mass model, stochastic disturbances without steps.
% 6.3: 6-state augmented model, stochastic disturbances with steps.
problem06Directory = fileparts(mfilename('fullpath'));
addpath(fullfile(problem06Directory,'..','common'),problem06Directory);

[problem06_2,problem06_2Figures] = problem06Simulate('6.2');
saveExperiment(fullfile(problem06Directory,'results'),'problem06_2', ...
    problem06_2,problem06_2Figures);

[problem06_3,problem06_3Figures] = problem06Simulate('6.3');
saveExperiment(fullfile(problem06Directory,'results'),'problem06_3', ...
    problem06_3,problem06_3Figures);
