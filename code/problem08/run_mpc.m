% Problem 8: unconstrained MPC on the linear four-tank model.
% Compares the 4-state (unbiased) and 6-state bias-augmented (biased) MPC
% on one reference step and two F3/F4 disturbance steps, and saves the
% metrics, MAT file, PNG and PDF figures in problem08/results/.
problem08Directory = fileparts(mfilename('fullpath'));
problem08CodeDirectory = fileparts(problem08Directory);
addpath(fullfile(problem08CodeDirectory,'common'), ...
    fullfile(problem08CodeDirectory,'problem04'), ...
    fullfile(problem08CodeDirectory,'problem06'), ...
    fullfile(problem08CodeDirectory,'problem07'), ...
    problem08Directory);

problem08Config = struct();
problem08Config.Ts = 5;
problem08Config.u0 = [300;300];
problem08Config.d0 = [20;20];
problem08Config.N = 20;
problem08Config.Wz = 10*eye(2);
problem08Config.Wu = 0.01*eye(2);
problem08Config.Wdu = 1*eye(2);
problem08Config.tFinal = 1000;
problem08Config.seed = 8001;
problem08Config.refTime = 200;
problem08Config.ref = [1;0.5];
problem08Config.stepTimes = [300 500];
problem08Config.stepSizes = [20;15];
problem08Config.settleTime = 100;
problem08Config.predictionTime = 700;

[problem08Result,problem08Figures] = problem08Simulate(problem08Config);
saveExperiment(fullfile(problem08Directory,'results'), ...
    'unconstrained_mpc',problem08Result,problem08Figures);

problem08Metrics = problem08Result.metrics;
problem08Names = {'unbiased','biased'};
fprintf(['\nProblem 8 metrics: tracking after t = %g s, final offset over ' ...
    'the last %g s\n'],problem08Config.refTime,problem08Config.settleTime);
for problem08Index = 1:numel(problem08Names)
    problem08Name = problem08Names{problem08Index};
    problem08M = problem08Metrics.(problem08Name);
    fprintf(['  %-9s RMSE [%7.3f %7.3f] cm  offset [%+7.3f %+7.3f] cm  ' ...
        'max|e| %6.3f cm  |u|^2 %8.1f  |du|^2 %8.1f\n'], ...
        problem08Name,problem08M.levelRMSE(1),problem08M.levelRMSE(2), ...
        problem08M.finalOffset(1),problem08M.finalOffset(2), ...
        max(problem08M.levelMaxAbs),problem08M.controlEffortTotal, ...
        problem08M.inputMoveEffortTotal);
end
fprintf('Problem 8 artifacts saved in %s\n', ...
    fullfile(problem08Directory,'results'));
