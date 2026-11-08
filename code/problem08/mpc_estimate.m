function [x_hat, w_hat, innovation] = mpc_estimate(mpc, x_prev, w_prev, u_prev, y)
%MPC_ESTIMATE One prediction and measurement update of the MPC estimator.
%
%   [x_hat, w_hat, innovation] = mpc_estimate(mpc, x_prev, w_prev, u_prev, y)
%
% Uses xp = A*x_prev + B*u_prev + G*w_prev, e = y - C*xp, then
% x_hat = xp + Kfx*e and w_hat = w_prev + Kfw*e (persistent) or Kfw*e
% (impulse). For nw = 0 the disturbance has no state and w_hat is empty.

if isempty(mpc.estimator)
    error('mpc_estimate:NoEstimator', ...
        'mpc_constructors was called without an estimator.');
end
validateattributes(x_prev,{'numeric'}, ...
    {'real','finite','size',[mpc.nx 1]},mfilename,'x_prev');
validateattributes(u_prev,{'numeric'}, ...
    {'real','finite','size',[mpc.nu 1]},mfilename,'u_prev');
validateattributes(y,{'numeric'}, ...
    {'real','finite','size',[mpc.nz 1]},mfilename,'y');

if mpc.nw == 0
    if ~isempty(w_prev)
        error('mpc_estimate:Dimension','w_prev must be empty when nw = 0.');
    end
    xp = mpc.A*x_prev + mpc.B*u_prev;
else
    validateattributes(w_prev,{'numeric'}, ...
        {'real','finite','size',[mpc.nw 1]},mfilename,'w_prev');
    xp = mpc.A*x_prev + mpc.B*u_prev + mpc.G*w_prev;
end

innovation = y - mpc.C*xp;
x_hat = xp + mpc.estimator.Kfx*innovation;
if mpc.nw == 0
    w_hat = [];
elseif mpc.estimator.persistent
    w_hat = w_prev + mpc.estimator.Kfw*innovation;
else
    w_hat = mpc.estimator.Kfw*innovation;
end
end
