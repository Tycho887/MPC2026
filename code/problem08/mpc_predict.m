function prediction = mpc_predict(mpc, x_hat, w_hat, U)
%MPC_PREDICT Roll out the prediction model over the horizon.
%
%   prediction = mpc_predict(mpc, x_hat, w_hat, U)
%
% Returns X (nx-by-(N+1)), Z (nz-by-N), U (nu-by-N), Zfree, b and Gamma.
% The stacked prediction is Z(:) = b + Gamma*U(:) with b = Phi_x*x + Phi_w*w.
% A persistent disturbance enters the recursion at every step; an impulse
% disturbance enters only the first step.

if nargin < 4 || isempty(U)
    error('mpc_predict:MissingInput','U is required.');
end
validateattributes(x_hat,{'numeric'},{'real','finite','size',[mpc.nx 1]}, ...
    mfilename,'x_hat');
if mpc.nw == 0
    if ~isempty(w_hat)
        error('mpc_predict:Dimension','w_hat must be empty when nw = 0.');
    end
else
    validateattributes(w_hat,{'numeric'}, ...
        {'real','finite','size',[mpc.nw 1]},mfilename,'w_hat');
end

Ucol = localStack(U,mpc.nu,mpc.N,'U');
Umat = reshape(Ucol,mpc.nu,mpc.N);

X = zeros(mpc.nx,mpc.N+1);
X(:,1) = x_hat;
Z = zeros(mpc.nz,mpc.N);
for j = 1:mpc.N
    xNext = mpc.A*X(:,j) + mpc.B*Umat(:,j);
    if mpc.nw > 0 && (mpc.persistent || j == 1)
        xNext = xNext + mpc.G*w_hat;
    end
    X(:,j+1) = xNext;
    Z(:,j) = mpc.C*xNext;
end

if mpc.nw == 0
    b = mpc.Phi_x*x_hat;
else
    b = mpc.Phi_x*x_hat + mpc.Phi_w*w_hat;
end
Zfree = reshape(b,mpc.nz,mpc.N);

Zstack = Z(:);
residual = Zstack - (b + mpc.Gamma*Ucol);
scale = max(1,norm(Zstack,inf));
if norm(residual,inf) > 1e-9*scale
    error('mpc_predict:InconsistentPrediction', ...
        'The recursion and the stacked prediction differ by %.3e.', ...
        norm(residual,inf));
end

prediction = struct();
prediction.X = X;
prediction.Z = Z;
prediction.U = Umat;
prediction.Zfree = Zfree;
prediction.b = b;
prediction.Gamma = mpc.Gamma;
end

function v = localStack(x,rows,cols,name)
validateattributes(x,{'numeric'},{'real','finite'},mfilename,name);
if isvector(x) && numel(x) == rows*cols
    v = x(:);
elseif isequal(size(x),[rows cols])
    v = x(:);
else
    error('mpc_predict:Dimension', ...
        '%s must have %d elements or size %d-by-%d.',name,rows*cols,rows,cols);
end
end
