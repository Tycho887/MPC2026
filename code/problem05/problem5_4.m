% Problems 5.4-5.5: exact ZOH models and pump-to-level Markov coefficients.
% Disturbance deviations are held constant over each sampling interval.
problem05Directory = fileparts(mfilename('fullpath'));
addpath(fullfile(problem05Directory,'..','common'));

%% Shared continuous-time model at the Problem 4 operating point
problem05Results = struct();
problem05Results.p = fourTankParameters();
problem05Results.Ts = 5;                % Sampling time [s]
problem05Results.count = 60;            % H_0 through H_59
problem05Results.op = fourTankSteadyState([300;300],[20;20],problem05Results.p);
problem05Results.continuous = linearizeFourTank(problem05Results.op,problem05Results.p);
% Measured and controlled outputs are both lower-tank levels.
problem05Results.continuous.Cz = problem05Results.continuous.C;

%% Two implementations of the same exact ZOH discretization
for problem05Method = {'zoh','expm'}
    problem05Name = problem05Method{1};
    problem05Model = discretizeLinearModel(problem05Results.continuous, ...
        problem05Results.Ts,problem05Name);
    % The fourth ss argument is zero feedthrough D, not the output matrix Cz.
    problem05System = ss(problem05Model.A,problem05Model.B, ...
        problem05Model.C,problem05Model.D,problem05Model.Ts);
    problem05Results.(problem05Name).model = problem05Model;
    problem05Results.(problem05Name).markov = markovParameters( ...
        problem05System,problem05Results.Ts,problem05Results.count);
end

problem05Results.differences = struct();
for problem05Field = {'A','B','E'}
    problem05Name = problem05Field{1};
    problem05Results.differences.(problem05Name) = max(abs( ...
        problem05Results.zoh.model.(problem05Name) - ...
        problem05Results.expm.model.(problem05Name)),[],'all');
end
problem05Results.differences.H = max(abs(problem05Results.zoh.markov.H - ...
    problem05Results.expm.markov.H),[],'all');
fprintf('ZOH/expm maximum absolute differences: A=%g, B=%g, E=%g, H=%g\n', ...
    problem05Results.differences.A,problem05Results.differences.B, ...
    problem05Results.differences.E,problem05Results.differences.H);

%% Markov coefficients: H_0 = D; H_k = C*Abar^(k-1)*Bbar for k >= 1
problem05Figures = experimentFigure('Name','Problem 5: ZOH Markov coefficients', ...
    'Position',[150 150 950 650]);
for problem05Output = 1:2
    for problem05Input = 1:2
        subplot(2,2,(problem05Output-1)*2+problem05Input);
        stem(problem05Results.zoh.markov.t, ...
            squeeze(problem05Results.zoh.markov.H(problem05Output,problem05Input,:)), ...
            'filled','MarkerSize',3,'Color',[0.1 0.4 0.8]);
        hold on;
        plot(problem05Results.expm.markov.t, ...
            squeeze(problem05Results.expm.markov.H(problem05Output,problem05Input,:)), ...
            '--','LineWidth',1.3,'Color',[0.85 0.35 0.1]);
        grid on;
        xlabel('Time kT_s (s)');
        ylabel(sprintf('H_{%d%d}(k)',problem05Output,problem05Input));
        title(sprintf('F_%d \\rightarrow h_%d (T_s = %g s)', ...
            problem05Input,problem05Output,problem05Results.Ts));
        legend('c2d ZOH','Matrix exponential','Location','best');
    end
end
saveExperiment(fullfile(problem05Directory,'results'),'discretization', ...
    problem05Results,problem05Figures);
