function [result, figures] = problem03Analysis()
% Migrated experiment; shared physical models live in common/.
figures = gobjects(0);
p = fourTankParameters();
u0 = [100;300]; d0 = [100;100];
op = fourTankSteadyState(u0,d0,p);
model = linearizeFourTank(op,p);
[Gm,Gh] = fourTankTransferFunctions(model);
% Step Response Simulation
% -------------------------------------------------------------------------
t_final = 1200; % Simulation time [s]
t_vec = linspace(0, t_final, 1000);

% Simulate unit step response for mass dynamics
[y_m, t_out_m] = step(Gm, t_vec);

% Simulate unit step response for height dynamics
[y_h, t_out_h] = step(Gh, t_vec);

% -------------------------------------------------------------------------
% Plotting Step Responses
% -------------------------------------------------------------------------
figures(end+1) = experimentFigure('Name', 'Four-Tank Step Responses', 'Color', 'w');

% Step in F1 -> [m1, m2]
subplot(2, 2, 1);
plot(t_out_m, y_m(:, 1, 1), 'b', 'LineWidth', 1.5);
grid on;
title('\delta F_1 \rightarrow \delta m_1 (Direct Path)');
xlabel('Time (s)'); ylabel('Mass change (g / (cm^3/s))');

subplot(2, 2, 3);
plot(t_out_m, y_m(:, 2, 1), 'r', 'LineWidth', 1.5);
grid on;
title('\delta F_1 \rightarrow \delta m_2 (Cross-coupled via Tank 4)');
xlabel('Time (s)'); ylabel('Mass change (g / (cm^3/s))');

% Step in F2 -> [m1, m2]
subplot(2, 2, 2);
plot(t_out_m, y_m(:, 1, 2), 'b', 'LineWidth', 1.5);
grid on;
title('\delta F_2 \rightarrow \delta m_1 (Cross-coupled via Tank 3)');
xlabel('Time (s)'); ylabel('Mass change (g / (cm^3/s))');

subplot(2, 2, 4);
plot(t_out_m, y_m(:, 2, 2), 'r', 'LineWidth', 1.5);
grid on;
title('\delta F_2 \rightarrow \delta m_2 (Direct Path)');
xlabel('Time (s)'); ylabel('Mass change (g / (cm^3/s))');

% -------------------------------------------------------------------------
% Static Decoupling and PI Controller Design
% -------------------------------------------------------------------------

controller = tuneFourTankPI(Gm);
T_closed = controller.closedLoop;
% -------------------------------------------------------------------------
% Plotting Closed-Loop Decoupled Step Responses
% -------------------------------------------------------------------------
figures(end+1) = experimentFigure('Name', 'Closed-Loop Decoupled Step Responses', 'Color', 'w');
step(T_closed, 1500);
title('Closed-Loop Step Response with Static Decoupling and PI Control');


result = struct('p',p,'op',op,'model',model,'Gm',Gm,'Gh',Gh, ...
    'controller',controller,'t',t_vec(:),'massSteps',y_m,'levelSteps',y_h);
end
