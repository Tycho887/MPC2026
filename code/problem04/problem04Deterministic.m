function [result, figures] = problem04Deterministic()
% Migrated experiment; shared physical models live in common/.
figures = gobjects(0);
p = fourTankParameters();
u0 = [300;300]; d = [20;20];
op = fourTankSteadyState(u0,d,p); x0 = op.x;
experiments = cell(2,3);
%% Simulation setup
t0 = 0;
t_step = 200;                    % Step time [s]
t_final = 1200;                  % End time [s]
step_fractions = [0.10, 0.25, 0.50];
line_colors = lines(numel(step_fractions));

% Run a separate set of step tests for each manipulated pump.
for pump_index = 1:2
    figures(end+1) = experimentFigure('Name', sprintf('Problem 4.1: Step response of F_%d', pump_index), ...
           'Color', 'w');

    ax_h1 = subplot(3,1,1);
    hold(ax_h1, 'on'); grid(ax_h1, 'on');
    ylabel(ax_h1, 'Level h_1 (cm)');
    title(ax_h1, sprintf('Tank levels after a step in pump F_%d', pump_index));

    ax_h2 = subplot(3,1,2);
    hold(ax_h2, 'on'); grid(ax_h2, 'on');
    ylabel(ax_h2, 'Level h_2 (cm)');

    ax_input = subplot(3,1,3);
    hold(ax_input, 'on'); grid(ax_input, 'on');
    ylabel(ax_input, sprintf('F_%d (cm^3/s)', pump_index));
    xlabel(ax_input, 'Time (s)');

    for step_index = 1:numel(step_fractions)
        fraction = step_fractions(step_index);
        u_step = u0;
        u_step(pump_index) = u0(pump_index) * (1 + fraction);

        % Integrate before and after the step with constant inputs on each
        % interval. Starting at the calculated steady state gives a clean,
        % deterministic step response without a startup transient.
        [t_pre, X_pre] = integrateFourTankInterval([t0,t_step],x0,u0,d,p);
        [t_post, X_post] = integrateFourTankInterval([t_step,t_final],X_pre(:,end),u_step,d,p);

        t_response = [t_pre; t_post(2:end)];
        X_response = [X_pre, X_post(:,2:end)];
        h_response = fourTankOutput(X_response,p)';
        input = repmat(u0,1,numel(t_response));
        input(pump_index,t_response>=t_step) = u_step(pump_index);
        experiments{pump_index,step_index} = struct('t',t_response,'x',X_response, ...
            'h',fourTankLevels(X_response,p),'z',h_response','y',h_response', ...
            'u',input,'d',repmat(d,1,numel(t_response)),'fraction',fraction);

        response_label = sprintf('+%.0f%%', 100 * fraction);
        input_label = sprintf('+%.0f%% (%.0f cm^3/s)', ...
            100 * fraction, u_step(pump_index));

        plot(ax_h1, t_response, h_response(:,1), ...
            'Color', line_colors(step_index,:), 'LineWidth', 1.5, ...
            'DisplayName', response_label);
        plot(ax_h2, t_response, h_response(:,2), ...
            'Color', line_colors(step_index,:), 'LineWidth', 1.5, ...
            'DisplayName', response_label);
        stairs(ax_input, [t0, t_step, t_final], ...
            [u0(pump_index), u_step(pump_index), u_step(pump_index)], ...
            'Color', line_colors(step_index,:), 'LineWidth', 1.5, ...
            'DisplayName', input_label);
    end

    xlim(ax_h1, [t0, t_final]);
    xlim(ax_h2, [t0, t_final]);
    xlim(ax_input, [t0, t_final]);
    legend(ax_h1, 'Location', 'best');
    legend(ax_h2, 'Location', 'best');
    legend(ax_input, 'Location', 'best');
end


result = struct('p',p,'op',op,'stepTime',t_step,'finalTime',t_final, ...
    'fractions',step_fractions,'solver','ode15s','experiments',{experiments});
end
