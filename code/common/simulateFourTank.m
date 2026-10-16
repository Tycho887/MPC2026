function result = simulateFourTank(t, x0, u, d, p, options)
%SIMULATEFOURTANK Open-loop ODE or additive SDE simulation.
% t: N-by-1 increasing grid; u,d: 2-by-N (held over each interval).
% options.solver: ode15s (default), ode45, or sde.
% SDE options: sigma (4-by-2) and dW (2-by-(N-1), including sqrt(dt)).
% options.measurementNoise: 2-by-N [cm], default zero.
% options.seed is metadata; callers generate random arrays explicitly.
if nargin < 6, options = struct(); end
if ~isfield(options,'solver'), options.solver = 'ode15s'; end
if ~isfield(options,'measurementNoise'), options.measurementNoise = zeros(2,numel(t)); end
t = t(:); N = numel(t);
assert(N >= 2 && all(diff(t)>0), 'Time grid must be strictly increasing.');
assert(isequal(size(u),[2 N]) && isequal(size(d),[2 N]), 'u and d must be 2-by-N.');
assert(isequal(size(options.measurementNoise),[2 N]), 'Noise must be 2-by-N.');
if strcmp(options.solver,'sde')
    assert(isequal(size(options.sigma),[4 2]) && ...
        isequal(size(options.dW),[2 N-1]), 'Invalid SDE noise dimensions.');
end
x = zeros(4,N); x(:,1) = x0;
for k = 1:N-1
    if strcmp(options.solver,'sde')
        % Constant diffusion makes the Milstein correction zero.
        dx = fourTankDynamics(t(k),x(:,k),u(:,k),d(:,k),p);
        x(:,k+1) = max(x(:,k) + diff(t(k:k+1))*dx + options.sigma*options.dW(:,k), 0);
    else
        [~, interval] = integrateFourTankInterval(t(k:k+1),x(:,k),u(:,k),d(:,k),p,options.solver);
        x(:,k+1) = interval(:,end);
    end
end
result = struct('t',t,'x',x,'h',fourTankLevels(x,p), ...
    'y',fourTankSensor(x,p,options.measurementNoise), ...
    'z',fourTankOutput(x,p),'u',u,'d',d,'p',p,'options',options);
end
