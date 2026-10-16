function [result, figures] = problem03LevelControl()
% Migrated experiment; shared physical models live in common/.
figures = gobjects(0);
p = fourTankParameters();
u0 = [300;300]; d0 = [100;100];
op = fourTankSteadyState(u0,d0,p); m_s = op.x;
model = linearizeFourTank(op,p);
[~,Gh] = fourTankTransferFunctions(model);
controller = tuneFourTankPI(Gh);
C1 = controller.C1; C2 = controller.C2; W = controller.W;
figures(end+1) = experimentFigure('Name','Linear level closed loop','Color','w');
step(controller.closedLoop,1500);
title('Linear level response with static decoupling and PI control');
%% 5. Nonlinear Simulation with Decoupled Discrete PI Control
Ts = 2.0;               % Sampling time [s]
t_final = 1200;         % Final simulation time [s]
t_vec = 0:Ts:t_final;
num_steps = length(t_vec);

h_target = [45; 30];    % Target operating point levels [cm]

% Discretize continuous controllers using Tustin (bilinear) approximation
C1_d = c2d(C1, Ts, 'tustin');
C2_d = c2d(C2, Ts, 'tustin');

% Convert to state-space for efficient discrete stepping inside the loop
sys1 = ss(C1_d);
sys2 = ss(C2_d);

xc1 = zeros(size(sys1.A, 1), 1);
xc2 = zeros(size(sys2.A, 1), 1);

x_current = m_s;        % Initialize states with steady-state masses
x_history = zeros(4,num_steps);
y_history = zeros(num_steps, 2);
u_history = zeros(num_steps, 2);

for k = 1:num_steps
    % Measure outputs (liquid levels in Tank 1 and Tank 2)
    h_current = fourTankOutput(x_current,p);
    x_history(:,k) = x_current;
    y_history(k, :) = fourTankSensor(x_current,p)';
    
    % Compute tracking error
    e_h = h_target - h_current;
    
    % Evaluate discrete controllers efficiently via state-space updates
    du1 = sys1.C * xc1 + sys1.D * e_h(1);
    du2 = sys2.C * xc2 + sys2.D * e_h(2);
    
    xc1 = sys1.A * xc1 + sys1.B * e_h(1);
    xc2 = sys2.A * xc2 + sys2.B * e_h(2);
    
    u_tilde = [du1; du2];
    
    % Apply static decoupling and add nominal steady-state flows
    u_physical = W * u_tilde;
    F_actual = u0 + u_physical;
    F_actual = max(F_actual, 0);     % Enforce non-negative actuator bounds
    u_history(k, :) = F_actual';
    
    % Integrate nonlinear plant over the sampling interval [t_k, t_{k+1}]
    if k < num_steps
        [~, X_ode] = integrateFourTankInterval(t_vec(k:k+1), x_current, F_actual, d0, p, 'ode45');
        x_current = X_ode(:,end);
    end
end

%% 6. Plot Nonlinear Simulation Results
figures(end+1) = experimentFigure('Name', 'Nonlinear Closed-Loop Simulation', 'Color', 'w');
subplot(2,1,1);
plot(t_vec, y_history, 'LineWidth', 1.5);
grid on;
legend('h_1', 'h_2', 'Location', 'Best');
title('Nonlinear Tank Levels under Decoupled Discrete PI Control');
xlabel('Time (s)'); ylabel('Level (cm)');

subplot(2,1,2);
plot(t_vec, u_history, 'LineWidth', 1.5);
grid on;
legend('F_1', 'F_2', 'Location', 'Best');
title('Manipulated Pump Flows');
xlabel('Time (s)'); ylabel('Flow (cm^3/s)');


result = struct('p',p,'op',op,'model',model,'controller',controller, ...
    't',t_vec(:),'x',x_history,'h',fourTankLevels(x_history,p), ...
    'y',y_history','z',fourTankOutput(x_history,p),'u',u_history', ...
    'd',repmat(d0,1,num_steps),'reference',repmat(h_target,1,num_steps), ...
    'sampleTime',Ts,'solver','ode45','pumpBounds',[0;Inf]);
end
