function [result, figures] = problem04Stochastic()
% Migrated experiment; shared physical models live in common/.
figures = gobjects(0);
p = fourTankParameters();
u0 = [300;300]; d_mean = [20;20];
op = fourTankSteadyState(u0,d_mean,p); x0 = op.x;
experiments = cell(2,3,3);
%% Stochastic simulation setup
dt = 0.1;                        % Integration step [s]
t0 = 0;
t_step = 200;                    % Step time [s]
t_final = 1200;                  % End time [s]
t = (t0:dt:t_final)';
N = numel(t);
step_fractions = [0.10, 0.25, 0.50];
line_colors = lines(numel(step_fractions));
noise_scales = [0.5, 1.0, 2.0];
noise_labels = {'Low (0.5x)', 'Medium (1x)', 'High (2x)'};

% Additive process noise enters the upper tanks, as in problem2_4_3.m.
sigma_F_base = [5; 5];           % Baseline diffusion [cm^3/s^(1/2)]


% Sensor noise is added to the measured levels of tanks 1 and 2.
sigma_y_base = 0.5;              % Baseline sensor standard deviation [cm]
R = diag([sigma_y_base^2, sigma_y_base^2]);
Lr = chol(R, 'lower');

% Reuse standardized noise paths across tests for fair comparisons.
rng(100, 'twister');
dW_unit = sqrt(dt) * randn(2, N - 1);
measurement_noise_unit = Lr * randn(2, N);
measurement_stride = 10;         % Plot sensor samples every 1 second
measurement_indices = 1:measurement_stride:N;

%% Simulate and plot a separate experiment for each pump
for pump_index = 1:2
    figures(end+1) = experimentFigure('Name', sprintf('Problem 4.2: F_%d step response', pump_index), ...
           'Units', 'inches', 'Position', [1, 1, 9, 9], 'Color', 'w');
    sgtitle(sprintf('Problem 4.2: Step responses to pump F_%d', pump_index));

    for noise_index = 1:numel(noise_scales)
        noise_scale = noise_scales(noise_index);
        sigma_mat = fourTankDiffusion(p,sigma_F_base * noise_scale);
        measurement_noise = noise_scale * measurement_noise_unit;

        ax_h1 = subplot(3,3,(noise_index - 1) * 3 + 1);
        hold(ax_h1, 'on'); grid(ax_h1, 'on');
        ylabel(ax_h1, 'h_1 (cm)');
        title(ax_h1, sprintf('%s noise: Tank 1', noise_labels{noise_index}));

        ax_h2 = subplot(3,3,(noise_index - 1) * 3 + 2);
        hold(ax_h2, 'on'); grid(ax_h2, 'on');
        ylabel(ax_h2, 'h_2 (cm)');
        title(ax_h2, sprintf('%s noise: Tank 2', noise_labels{noise_index}));

        ax_input = subplot(3,3,(noise_index - 1) * 3 + 3);
        hold(ax_input, 'on'); grid(ax_input, 'on');
        ylabel(ax_input, sprintf('F_%d (cm^3/s)', pump_index));
        title(ax_input, 'Pump input');

        for step_index = 1:numel(step_fractions)
            fraction = step_fractions(step_index);
            u_step = u0;
            u_step(pump_index) = u0(pump_index) * (1 + fraction);

            options = struct('solver','sde','sigma',sigma_mat,'dW',dW_unit, ...
                'measurementNoise',measurement_noise,'seed',100);
            u_profile = repmat(u0,1,N); u_profile(:,t>=t_step) = repmat(u_step,1,sum(t>=t_step));
            r = simulateFourTank(t,x0,u_profile,repmat(d_mean,1,N),p,options);
            h_true = r.z; h_measured = r.y; u_history = r.u;
            experiments{pump_index,noise_index,step_index} = r;

            response_label = sprintf('+%.0f%%', 100 * fraction);

            plot(ax_h1, t, h_true(1,:), ...
                'Color', line_colors(step_index,:), 'LineWidth', 1.2, ...
                'DisplayName', response_label);
            plot(ax_h1, t(measurement_indices), ...
                h_measured(1,measurement_indices), '.', ...
                'Color', line_colors(step_index,:), 'MarkerSize', 4, ...
                'HandleVisibility', 'off');

            plot(ax_h2, t, h_true(2,:), ...
                'Color', line_colors(step_index,:), 'LineWidth', 1.2, ...
                'DisplayName', response_label);
            plot(ax_h2, t(measurement_indices), ...
                h_measured(2,measurement_indices), '.', ...
                'Color', line_colors(step_index,:), 'MarkerSize', 4, ...
                'HandleVisibility', 'off');

            stairs(ax_input, t, u_history(pump_index,:), ...
                'Color', line_colors(step_index,:), 'LineWidth', 1.2, ...
                'DisplayName', response_label);
        end

        xlim(ax_h1, [t0, t_final]);
        xlim(ax_h2, [t0, t_final]);
        xlim(ax_input, [t0, t_final]);
        xlabel(ax_h1, 'Time (s)');
        xlabel(ax_h2, 'Time (s)');
        xlabel(ax_input, 'Time (s)');

        if noise_index == 1
            % Step colors identify the 10%, 25%, and 50% input changes.
            % Measurement markers use the same color as their true response.
            plot(ax_h1, NaN, NaN, 'k.', 'MarkerSize', 9, ...
                'DisplayName', 'Measurements (1 s shown)');
            legend(ax_h1, 'Location', 'best');
            legend(ax_input, 'Location', 'best');
        end
    end
end


result = struct('p',p,'op',op,'stepTime',t_step,'finalTime',t_final, ...
    'fractions',step_fractions,'noiseScales',noise_scales,'seed',100,'experiments',{experiments});
end
