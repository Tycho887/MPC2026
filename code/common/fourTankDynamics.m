function dx = fourTankDynamics(~, x, u, d, p)
%FOURTANKDYNAMICS Mass balance dx/dt=f(t,x,u,d,p) [g/s].
% The same drift serves deterministic, held-disturbance, and SDE models.
[Q, U, D] = fourTankRouting(p);
levels = fourTankLevels(max(x, 0), p); % Protect solver trial evaluations.
q = p.a .* sqrt(2 * p.g * levels);
dx = p.rho * (Q * q + U * u + D * d);
end
