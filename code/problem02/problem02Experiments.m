function [results, figures] = problem02Experiments()
p = fourTankParameters();
results = struct(); figures = gobjects(0);
names = {'deterministic','piecewise_stochastic','sde'};
for mode = 1:3
    if mode == 3, dt = 0.1; else, dt = 10; end
    t = (0:dt:1200)'; N = numel(t);
    u = 300*ones(2,N);
    u(1,floor(N/3):end) = 350;
    u(2,floor(2*N/3):end) = 250;
    d = repmat([20;20],1,N);
    options = struct('solver','ode15s','seed',100);
    rng(options.seed,'twister');
    if mode == 2
        d = d + diag([5;5])*randn(2,N);
        options.measurementNoise = 0.5*randn(2,N);
    elseif mode == 3
        options.solver = 'sde';
        options.sigma = fourTankDiffusion(p,[5;5]);
        options.dW = sqrt(dt)*randn(2,N-1);
        options.measurementNoise = 0.5*randn(2,N);
    end
    r = simulateFourTank(t,zeros(4,1),u,d,p,options);
    results.(names{mode}) = r;
    figures(end+1) = experimentFigure('Name',names{mode},'Color','w');
    subplot(3,1,1); stairs(t,r.u','LineWidth',1.5); grid on;
    ylabel('Pump flow (cm^3/s)'); legend('F_1','F_2'); title(strrep(names{mode},'_',' '));
    subplot(3,1,2);
    if mode == 3
        plot(t,r.h(3:4,:)'); ylabel('Level (cm)'); legend('h_3','h_4');
    else
        stairs(t,r.d'); ylabel('Disturbance (cm^3/s)'); legend('F_3','F_4');
    end
    grid on;
    subplot(3,1,3); plot(t,r.y','.','MarkerSize',3); hold on;
    plot(t,r.z','LineWidth',1.2); grid on;
    ylabel('Level (cm)'); xlabel('Time (s)'); legend('y_1','y_2','z_1','z_2');
end
end
