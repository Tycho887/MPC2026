function [t, x] = integrateFourTankInterval(tspan, x0, u, d, p, solver)
%INTEGRATEFOURTANKINTERVAL Integrate with held inputs; x is 4-by-numel(t).
if nargin < 6, solver = 'ode15s'; end
assert(ismember(solver, {'ode15s', 'ode45'}), 'Unsupported ODE solver.');
[t, states] = feval(solver, @(time,state) fourTankDynamics(time,state,u,d,p), ...
    tspan, x0);
x = states';
end
