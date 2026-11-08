function [x_ss, u_ss] = mpc_target(mpc, r, w_hat)
%MPC_TARGET Steady-state state and input target for a constant reference.
%
%   [x_ss, u_ss] = mpc_target(mpc, r, w_hat)
%
% Solves M*u_ss = r - C*(I-A)\G*w with M = C*(I-A)\B and returns
% x_ss = (I-A)\(B*u_ss + G*w). An empty w_hat is treated as a zero
% disturbance when nw > 0.

validateattributes(r,{'numeric'},{'real','finite'},mfilename,'r');
if numel(r) ~= mpc.nz
    error('mpc_target:Dimension','r must have %d elements.',mpc.nz);
end
r = r(:);

if mpc.nw == 0
    if ~isempty(w_hat)
        error('mpc_target:Dimension','w_hat must be empty when nw = 0.');
    end
    disturbance = zeros(mpc.nx,1);
else
    if isempty(w_hat)
        w = zeros(mpc.nw,1);
    else
        validateattributes(w_hat,{'numeric'}, ...
            {'real','finite','size',[mpc.nw 1]},mfilename,'w_hat');
        w = w_hat;
    end
    if mpc.persistent
        disturbance = mpc.G*w;
    else
        disturbance = zeros(mpc.nx,1);
    end
end

IminusA = eye(mpc.nx) - mpc.A;
if rcond(IminusA) < 1e-12
    error('mpc_target:SingularDynamics', ...
        'I - A is singular; a steady-state target does not exist.');
end
resolventB = IminusA\mpc.B;
M = mpc.C*resolventB;
if rank(M) < mpc.nu
    error('mpc_target:UnreachableReference', ...
        'M = C*(I-A)\B does not have full column rank.');
end

u_ss = M\(r - mpc.C*(IminusA\disturbance));
x_ss = IminusA\(mpc.B*u_ss + disturbance);

outputResidual = mpc.C*x_ss - r;
if norm(outputResidual,inf) > 1e-8*max(1,norm(r,inf))
    warning('mpc_target:ReferenceResidual', ...
        'The target output misses the reference by %.3e.', ...
        norm(outputResidual,inf));
end
stateResidual = x_ss - (mpc.A*x_ss + mpc.B*u_ss + disturbance);
if norm(stateResidual,inf) > 1e-8*max(1,norm(x_ss,inf))
    warning('mpc_target:StateResidual', ...
        'The target state misses the steady-state equation by %.3e.', ...
        norm(stateResidual,inf));
end
end
