function [u_k, x_hat, w_hat, info] = mpc_step(mpc, x_hat_prev, w_hat_prev, ...
    u_prev, y_k, R_k, U_bar_k)
%MPC_STEP One unconstrained MPC step: estimate, predict and solve.
%
%   [u_k, x_hat, w_hat, info] = mpc_step(mpc, x_hat_prev, w_hat_prev, ...
%       u_prev, y_k, R_k, U_bar_k)
%
% R_k is the stacked reference over the horizon (nz*N elements or nz-by-N).
% U_bar_k is the stacked steady-state input target (nu*N elements or
% nu-by-N). When U_bar_k is empty it is computed from mpc_target using the
% first reference block and the updated disturbance estimate.
%
% info contains U, Z, b, c, g, rho, innovation, prediction and the QP data.
% The online problem is solved with qpsolver from Problem 7. Problem 8 passes
% infinite variable bounds and an empty constraint matrix; later problems can
% supply finite bounds through mpc.l_bar/mpc.u_bar and mpc.A_bar/mpc.bl_bar/
% mpc.bu_bar without changing this function.

if nargin < 7
    U_bar_k = [];
end
validateattributes(x_hat_prev,{'numeric'}, ...
    {'real','finite','size',[mpc.nx 1]},mfilename,'x_hat_prev');
validateattributes(u_prev,{'numeric'}, ...
    {'real','finite','size',[mpc.nu 1]},mfilename,'u_prev');
validateattributes(y_k,{'numeric'}, ...
    {'real','finite','size',[mpc.nz 1]},mfilename,'y_k');
if mpc.nw == 0
    if ~isempty(w_hat_prev)
        error('mpc_step:Dimension','w_hat_prev must be empty when nw = 0.');
    end
else
    validateattributes(w_hat_prev,{'numeric'}, ...
        {'real','finite','size',[mpc.nw 1]},mfilename,'w_hat_prev');
end

[x_hat, w_hat, innovation] = mpc_estimate(mpc, x_hat_prev, w_hat_prev, ...
    u_prev, y_k);

R_k = localStack(R_k,mpc.nz,mpc.N,'R_k');
if mpc.nw == 0
    b = mpc.Phi_x*x_hat;
else
    b = mpc.Phi_x*x_hat + mpc.Phi_w*w_hat;
end
c = R_k - b;

if isempty(U_bar_k)
    r = R_k(1:mpc.nz);
    [~, u_ss] = mpc_target(mpc, r, w_hat);
    U_bar = repmat(u_ss, mpc.N, 1);
else
    U_bar = localStack(U_bar_k,mpc.nu,mpc.N,'U_bar_k');
end

g = mpc.Mz*c + mpc.Mu*U_bar + mpc.Mdu*u_prev;
if exist('qpsolver','file') ~= 2
    error('mpc_step:MissingQPSolver', ...
        'qpsolver is not on the path; add code/problem07.');
end
[U,solverInfo] = qpsolver(mpc.H,g,mpc.l_bar,mpc.u_bar, ...
    mpc.A_bar,mpc.bl_bar,mpc.bu_bar,mpc.xinit);
if isempty(U) || any(~isfinite(U))
    error('mpc_step:InvalidSolution', ...
        'qpsolver returned no finite solution (exitflag %d).',solverInfo);
end
if solverInfo ~= 1
    warning('mpc_step:QPSolverExit', ...
        'qpsolver returned exitflag %d.',solverInfo);
end
residual = mpc.H*U + g;
if norm(residual,inf) > 1e-6*max(1,norm(g,inf))
    warning('mpc_step:SolveResidual', ...
        'H*U + g has infinity norm %.3e.',norm(residual,inf));
end

wc = mpc.Wz_bar*c;
wu = mpc.Wu_bar*U_bar;
wdu = mpc.Wdu_bar*(mpc.I0*u_prev);
rho = 0.5*(c'*wc + U_bar'*wu + u_prev'*(mpc.I0'*wdu));

u_k = U(1:mpc.nu);
prediction = mpc_predict(mpc,x_hat,w_hat,U);

info = struct();
info.U = U;
info.Z = prediction.Z;
info.b = b;
info.c = c;
info.g = g;
info.rho = rho;
info.innovation = innovation;
info.prediction = prediction;
info.solverInfo = solverInfo;
info.qp = struct('H',mpc.H,'g',g,'l',mpc.l_bar,'u',mpc.u_bar, ...
    'A',mpc.A_bar,'bl',mpc.bl_bar,'bu',mpc.bu_bar, ...
    'xinit',mpc.xinit,'exitflag',solverInfo);
end

function v = localStack(x,rows,cols,name)
validateattributes(x,{'numeric'},{'real','finite'},mfilename,name);
if isvector(x) && numel(x) == rows*cols
    v = x(:);
elseif isequal(size(x),[rows cols])
    v = x(:);
else
    error('mpc_step:Dimension', ...
        '%s must have %d elements or size %d-by-%d.',name,rows*cols,rows,cols);
end
end
