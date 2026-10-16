% Problem 2: deterministic, held stochastic disturbances, and additive SDE.
problemDirectory = fileparts(mfilename('fullpath'));
addpath(fullfile(problemDirectory,'..','common'), problemDirectory);
[problem02Results, problem02Figures] = problem02Experiments();
saveExperiment(fullfile(problemDirectory,'results'),'open_loop',problem02Results,problem02Figures);
