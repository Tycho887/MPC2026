function G = fitContinuousStepModel(t_rel, y_step, order)
    % Least-squares continuous-time step response fitting
    if order == 1
        % Model: y(t) = K * (1 - exp(-t / tau))
        theta0 = [mean(y_step(end-20:end)), 60]; % Initial guess [K, tau]
        obj = @(th) sum((y_step - th(1) * (1 - exp(-max(t_rel, 0) / max(th(2), 1e-2)))).^2);
        opt = optimset('Display', 'off');
        th_est = fminsearch(obj, theta0, opt);
        K = th_est(1);
        tau = max(th_est(2), 1e-2);
        G = tf(K, [tau, 1]);
    else
        % Model: y(t) = K * (1 - (tau1*exp(-t/tau1) - tau2*exp(-t/tau2))/(tau1 - tau2))
        theta0 = [mean(y_step(end-20:end)), 40, 80]; % Initial guess [K, tau1, tau2]
        obj = @(th) sum((y_step - evalSecondOrderStep(t_rel, th(1), th(2), th(3))).^2);
        opt = optimset('Display', 'off');
        th_est = fminsearch(obj, theta0, opt);
        K = th_est(1);
        tau1 = max(th_est(2), 1e-2);
        tau2 = max(th_est(3), 1e-2);
        G = tf(K, conv([tau1, 1], [tau2, 1]));
    end
end

function y = evalSecondOrderStep(t, K, t1, t2)
    t = max(t, 0);
    t1 = max(t1, 1e-2);
    t2 = max(t2, 1e-2);
    if abs(t1 - t2) < 1e-3
        t2 = t1 + 1e-3;
    end
    y = K * (1 - (t1 * exp(-t / t1) - t2 * exp(-t / t2)) / (t1 - t2));
end
