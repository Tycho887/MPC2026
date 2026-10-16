function [result, figures] = problem04Identification(backend)
% Migrated experiment; shared physical models live in common/.
figures = gobjects(0);
p = fourTankParameters();
if nargin < 1, backend = 'auto'; end
u0 = [300;300]; d_mean = [20;20];
op = fourTankSteadyState(u0,d_mean,p); x0 = op.x; h_s = op.h;
experiments = cell(1,2); normalized = cell(2,2);
dt = 0.1;                        % Simulation integration step [s]
t_step = 200;                    % Step time [s]
t_final = 1200;                  % End time [s]
t = (0:dt:t_final)';
N = numel(t);

step_fraction = 0.25;            % +25% step experiment
sigma_F = [5; 5];                % Process diffusion
sigma_mat = fourTankDiffusion(p,sigma_F);
sigma_y = 0.5;                   % Measurement noise std dev [cm]
Lr = diag([sigma_y, sigma_y]);

% Identification downsampling (1 second interval reduces noise overfit)
sample_decim = round(1.0 / dt);
t_id = t(1:sample_decim:end);

has_ident_toolbox = strcmp(backend,'toolbox') || (strcmp(backend,'auto') && ...
    exist('tfest','file') ~= 0 && exist('iddata','file') ~= 0 && license('test','Identification_Toolbox'));
assert(ismember(backend,{'auto','toolbox','least_squares'}),'Unknown identification backend.');

%% 3. Execute Step Experiments and Identify Channels
% Transfer function container: G(i,j) corresponds to input j -> output i
G_identified = cell(2,2);
model_orders = [1, 2; 2, 1];     % [G11, G12; G21, G22] physical pole counts

figures(end+1) = experimentFigure('Name', 'Problem 4.4: Normalized Step Responses and Fitted Models', ...
       'Color', 'w', 'Position', [100, 100, 1050, 750]);

for j_in = 1:2
    % Generate reproducible noise realization
    rng(100 + j_in, 'twister');
    dW = sqrt(dt) * randn(2, N - 1);
    measurement_noise = Lr * randn(2, N);

    % Configure step input
    delta_u = step_fraction * u0(j_in);
    u_step = u0;
    u_step(j_in) = u0(j_in) + delta_u;

    % Simulate stochastic plant
    options = struct('solver','sde','sigma',sigma_mat,'dW',dW, ...
        'measurementNoise',measurement_noise,'seed',100+j_in);
    u_profile = repmat(u0,1,N); u_profile(:,t>=t_step) = repmat(u_step,1,sum(t>=t_step));
    r = simulateFourTank(t,x0,u_profile,repmat(d_mean,1,N),p,options);
    experiments{j_in} = r; h_true = r.z; h_meas = r.y;

    % Baseline pre-step level estimation
    pre_idx = (t >= (t_step - 50)) & (t < t_step);
    h0_est = mean(h_meas(:, pre_idx), 2);

    for i_out = 1:2
        % Output deviation scaled by input step amplitude
        delta_h_meas = h_meas(i_out, :)' - h0_est(i_out);
        y_norm_meas = delta_h_meas / delta_u;

        delta_h_true = h_true(i_out, :)' - h_s(i_out);
        y_norm_true = delta_h_true / delta_u;

        % Subsampled data for identification
        y_norm_id = y_norm_meas(1:sample_decim:end);

        order_ij = model_orders(i_out, j_in);

        G_identified{i_out,j_in} = identifyFourTankChannel(t_id,y_norm_id,t_step,order_ij,has_ident_toolbox);
        normalized{i_out,j_in} = struct('t',t,'true',y_norm_true,'measured',y_norm_meas, ...
            'identificationTime',t_id,'identificationOutput',y_norm_id, ...
            'sampleTime',dt,'identificationSampleTime',t_id(2)-t_id(1),'baseline',h0_est(i_out),'inputStep',delta_u);

        % Model simulation for validation
        u_unit = double(t >= t_step);
        y_model = lsim(G_identified{i_out, j_in}, u_unit, t);

        normalized{i_out,j_in}.model = y_model;

        % Plotting
        subplot_idx = (i_out - 1) * 2 + j_in;
        subplot(2, 2, subplot_idx);
        plot(t, y_norm_true, 'k-', 'LineWidth', 1.2, 'DisplayName', 'True Nonlin. Step');
        hold on; grid on;
        plot(t_id, y_norm_id, 'b.', 'MarkerSize', 4, 'DisplayName', 'Noisy Meas. / \Delta u');
        plot(t, y_model, 'r--', 'LineWidth', 1.8, 'DisplayName', sprintf('Identified (%d-pole)', order_ij));

        xlabel('Time (s)');
        ylabel(sprintf('S_{%d%d} [cm / (cm^3/s)]', i_out, j_in));
        title(sprintf('Channel F_%d \\rightarrow h_%d (Order %d)', j_in, i_out, order_ij));
        if subplot_idx == 1
            legend('Location', 'southeast');
        end
    end
end

%% 4. Summary of Identified Models
fprintf('====================================================================\n');
fprintf('             IDENTIFIED CONTINUOUS-TIME TRANSFER FUNCTIONS          \n');
fprintf('====================================================================\n\n');
for i_out = 1:2
    for j_in = 1:2
        fprintf('--- Transfer Function G_{%d%d}(s) (F_%d -> h_%d) ---\n', ...
            i_out, j_in, j_in, i_out);
        sys = minreal(G_identified{i_out, j_in});
        [num, den] = tfdata(sys, 'v');
        fprintf('Numerator: %s\nDenominator: %s\n',mat2str(num,6),mat2str(den,6));
        gain_k = num(end) / den(end);
        poles_p = roots(den);
        time_constants = -1 ./ real(poles_p);
        fprintf('Static Gain K: %.4f cm / (cm^3/s)\n', gain_k);
        fprintf('Time Constant(s) tau: %s seconds\n\n', ...
            mat2str(round(time_constants', 2)));
    end
end


G = [G_identified{1,1} G_identified{1,2};G_identified{2,1} G_identified{2,2}];
if has_ident_toolbox, usedBackend = 'toolbox'; else, usedBackend = 'least_squares'; end
result = struct('p',p,'op',op,'G',G,'channels',{G_identified}, ...
    'modelOrders',model_orders,'backend',usedBackend,'stepTime',t_step, ...
    'fraction',step_fraction,'experiments',{experiments},'normalized',{normalized});
end
