function mpc = mpc_constructors(model, config, estimator)
%MPC_CONSTRUCTORS Build an unconstrained deviation-variable MPC.
%
%   mpc = mpc_constructors(model, config, estimator)
%
% model contains A, B, C and optional G, D and Ts. D must be zero. config
% contains N, Wz, Wu, Wdu and disturbanceMode ('persistent', 'impulse' or
% 'none'). The optional estimator contains the Kalman gains Kfx and Kfw, its
% disturbance dimension nw and the persistent flag; these set the prediction
% mode, while config.disturbanceMode supplies it when no estimator is given.
% Weights may be scalars, vectors or square matrices.
%
% The decision variable is the stacked input U = [u_k; ...; u_{k+N-1}]. The
% prediction is Z = Phi_x*x + Phi_w*w + Gamma*U, and the cost uses the
% steady-state input target U_bar. The online step solves the standard-form
% QP of the Problem 7 qpsolver interface, with g = Mz*c + Mu*U_bar + Mdu*u_{k-1}
% and c = R - b. Problem 8 supplies infinite variable bounds and no linear
% constraints; the same fields carry finite bounds in later problems.

if nargin < 1 || isempty(model)
    error('mpc_constructors:MissingModel','model is required.');
end
if nargin < 2 || isempty(config)
    config = struct();
end
if nargin < 3
    estimator = [];
end

A = localModelField(model,'A');
C = localModelField(model,'C');
[nx,nxCheck] = size(A);
if nxCheck ~= nx
    error('mpc_constructors:ModelDimension','model.A must be square.');
end
nz = size(C,1);
if size(C,2) ~= nx
    error('mpc_constructors:ModelDimension','model.C must have %d columns.',nx);
end

B = localModelField(model,'B');
nu = size(B,2);
if size(B,1) ~= nx
    error('mpc_constructors:ModelDimension','model.B must have %d rows.',nx);
end

D = zeros(nz,nu);
if isfield(model,'D') && ~isempty(model.D)
    D = model.D;
    validateattributes(D,{'numeric'},{'real','finite'},mfilename,'model.D');
    if isscalar(D) && D == 0
        D = zeros(nz,nu);
    else
        if ~isequal(size(D),[nz nu])
            error('mpc_constructors:ModelDimension', ...
                'model.D must be %d-by-%d.',nz,nu);
        end
        if any(D(:) ~= 0)
            error('mpc_constructors:NonzeroD','model.D must be zero.');
        end
    end
end

if isfield(model,'G') && ~isempty(model.G)
    G = localModelField(model,'G');
    if size(G,1) ~= nx
        error('mpc_constructors:ModelDimension','model.G must have %d rows.',nx);
    end
else
    G = zeros(nx,0);
end

Ts = [];
if isfield(model,'Ts') && ~isempty(model.Ts)
    validateattributes(model.Ts,{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,'model.Ts');
    Ts = model.Ts;
end

config = localDefault(config,'N',20);
config = localDefault(config,'Wz',eye(nz));
config = localDefault(config,'Wu',1e-2*eye(nu));
config = localDefault(config,'Wdu',1e-1*eye(nu));
validateattributes(config.N,{'numeric'}, ...
    {'real','finite','scalar','integer','positive'},mfilename,'config.N');
N = double(config.N);

Wz = localWeight(config.Wz,nz,'config.Wz');
Wu = localWeight(config.Wu,nu,'config.Wu');
Wdu = localWeight(config.Wdu,nu,'config.Wdu');

%% Prediction mode and estimator gains
hasEstimator = ~isempty(estimator);
nw = 0;
isPersistent = false;
if hasEstimator
    if ~isfield(estimator,'nw') || isempty(estimator.nw)
        error('mpc_constructors:EstimatorField','estimator.nw is required.');
    end
    validateattributes(estimator.nw,{'numeric'}, ...
        {'real','finite','scalar','integer','nonnegative'}, ...
        mfilename,'estimator.nw');
    nw = double(estimator.nw);
    if isfield(estimator,'persistent') && ~isempty(estimator.persistent)
        validateattributes(estimator.persistent,{'logical','numeric'}, ...
            {'real','scalar'},mfilename,'estimator.persistent');
        isPersistent = logical(estimator.persistent);
    end
    if nw == 0 && isPersistent
        error('mpc_constructors:ModeMismatch', ...
            'estimator.persistent requires estimator.nw > 0.');
    end
end

hasMode = isfield(config,'disturbanceMode') && ~isempty(config.disturbanceMode);
if hasMode
    mode = validatestring(config.disturbanceMode, ...
        {'persistent','impulse','none'},mfilename,'config.disturbanceMode');
else
    mode = '';
end

if hasEstimator
    if nw == 0
        if hasMode && ~strcmp(mode,'none')
            warning('mpc_constructors:ModeOverride', ...
                'disturbanceMode = %s ignored: estimator.nw = 0.',mode);
        end
        mode = 'none';
        isPersistent = false;
    else
        if size(G,2) ~= nw
            error('mpc_constructors:DisturbanceDimension', ...
                'estimator.nw = %d requires model.G with %d columns.',nw,nw);
        end
        if isPersistent
            estimatorMode = 'persistent';
        else
            estimatorMode = 'impulse';
        end
        if hasMode && ~strcmp(mode,estimatorMode)
            warning('mpc_constructors:ModeOverride', ...
                'disturbanceMode = %s overridden as %s by the estimator.', ...
                mode,estimatorMode);
        end
        mode = estimatorMode;
    end
else
    if ~hasMode || strcmp(mode,'none')
        nw = 0;
        isPersistent = false;
        mode = 'none';
    else
        if isempty(G) || size(G,2) == 0
            warning('mpc_constructors:ModeOverride', ...
                'disturbanceMode = %s ignored: model.G is empty.',mode);
            nw = 0;
            isPersistent = false;
            mode = 'none';
        else
            nw = size(G,2);
            isPersistent = strcmp(mode,'persistent');
        end
    end
end

est = [];
if hasEstimator
    if ~isfield(estimator,'Kfx') || isempty(estimator.Kfx)
        error('mpc_constructors:EstimatorField','estimator.Kfx is required.');
    end
    Kfx = estimator.Kfx;
    validateattributes(Kfx,{'numeric'},{'real','finite','size',[nx nz]}, ...
        mfilename,'estimator.Kfx');
    if nw > 0
        if ~isfield(estimator,'Kfw') || isempty(estimator.Kfw)
            error('mpc_constructors:EstimatorField','estimator.Kfw is required.');
        end
        Kfw = estimator.Kfw;
        validateattributes(Kfw,{'numeric'},{'real','finite','size',[nw nz]}, ...
            mfilename,'estimator.Kfw');
    else
        if isfield(estimator,'Kfw') && ~isempty(estimator.Kfw)
            error('mpc_constructors:EstimatorField', ...
                'estimator.Kfw must be empty when estimator.nw = 0.');
        end
        Kfw = [];
    end

    if nw > 0
        w0 = zeros(nw,1);
    else
        w0 = [];
    end
    est = estimator;
    est = localDefault(est,'kind','custom');
    est = localDefault(est,'mode','');
    est = localDefault(est,'x0',zeros(nx,1));
    est = localDefault(est,'w0',w0);
    est = localDefault(est,'model',model);
    est = localDefault(est,'filter',[]);
    est.nx = nx;
    est.nw = nw;
    est.Kfx = double(Kfx);
    est.Kfw = double(Kfw);
    est.persistent = isPersistent;
end

%% Stacked prediction matrices
Gamma = zeros(nz*N,nu*N);
Phi_x = zeros(nz*N,nx);
Phi_w = zeros(nz*N,nw);

powers = cell(N+1,1);
powers{1} = eye(nx);
for k = 1:N
    powers{k+1} = A*powers{k};
end

powerSum = zeros(nx);
for j = 1:N
    rows = (j-1)*nz + (1:nz);
    Phi_x(rows,:) = C*powers{j+1};
    for i = 1:j
        cols = (i-1)*nu + (1:nu);
        Gamma(rows,cols) = C*powers{j-i+1}*B;
    end
    if nw > 0
        if isPersistent
            powerSum = powerSum + powers{j};
            Phi_w(rows,:) = C*powerSum*G;
        else
            Phi_w(rows,:) = C*powers{j}*G;
        end
    end
end

%% Stacked weights and difference operator
Wz_bar = kron(eye(N),Wz);
Wu_bar = kron(eye(N),Wu);
Wdu_bar = kron(eye(N),Wdu);

I0 = [eye(nu); zeros(nu*(N-1),nu)];
Lambda = eye(nu*N) - diag(ones(nu*(N-1),1),-nu);

Wz_sq = Wz_bar'*Wz_bar;
Wu_sq = Wu_bar'*Wu_bar;
Wdu_sq = Wdu_bar'*Wdu_bar;
H = Gamma'*Wz_sq*Gamma + Wu_sq + Lambda'*Wdu_sq*Lambda;
H = (H + H')/2;
Mz = -Gamma'*Wz_sq;
Mu = -Wu_sq;
Mdu = -Lambda'*Wdu_sq*I0;

%% Validation and factorization
sv = svd(Gamma);
tolRank = max(size(Gamma))*eps(max(sv));
if numel(sv) < nu*N || sum(sv > tolRank) < nu*N
    error('mpc_constructors:RankDeficientGamma', ...
        'Gamma must have full column rank; increase N or check controllability.');
end
if any(~isfinite(H(:)))
    error('mpc_constructors:InvalidH','H must be finite.');
end

eigH = eig(H);
scale = max(1,max(abs(eigH)));
if min(eigH) < -1e-10*scale
    error('mpc_constructors:IndefiniteH','H must be positive semidefinite.');
end
Hpd = all(eigH > 0);
if ~Hpd
    warning('mpc_constructors:SingularH', ...
        'H is not positive definite; the QP may have multiple minimizers.');
end

mpc = struct();
mpc.A = A;
mpc.B = B;
mpc.C = C;
mpc.G = G;
mpc.D = D;
mpc.Ts = Ts;
mpc.nx = nx;
mpc.nu = nu;
mpc.nz = nz;
mpc.nw = nw;
mpc.N = N;
mpc.disturbanceMode = mode;
mpc.persistent = isPersistent;
mpc.Wz = Wz;
mpc.Wu = Wu;
mpc.Wdu = Wdu;
mpc.Wz_bar = Wz_bar;
mpc.Wu_bar = Wu_bar;
mpc.Wdu_bar = Wdu_bar;
mpc.Lambda = Lambda;
mpc.I0 = I0;
mpc.Gamma = Gamma;
mpc.Phi_x = Phi_x;
mpc.Phi_w = Phi_w;
mpc.H = H;
mpc.Mz = Mz;
mpc.Mu = Mu;
mpc.Mdu = Mdu;
mpc.Hpd = Hpd;
mpc.nU = nu*N;
mpc.l_bar = -inf(nu*N,1);
mpc.u_bar = inf(nu*N,1);
mpc.A_bar = zeros(0,nu*N);
mpc.bl_bar = zeros(0,1);
mpc.bu_bar = zeros(0,1);
mpc.xinit = zeros(nu*N,1);
mpc.estimator = est;
mpc.model = model;
mpc.config = config;
end

function M = localModelField(model,name)
if ~isfield(model,name) || isempty(model.(name))
    error('mpc_constructors:ModelField','model.%s is required.',name);
end
M = model.(name);
validateattributes(M,{'numeric'},{'real','finite'},mfilename,['model.' name]);
end

function W = localWeight(W,dim,name)
validateattributes(W,{'numeric'},{'real','finite'},mfilename,name);
if isscalar(W)
    W = double(W)*eye(dim);
elseif isvector(W)
    if numel(W) ~= dim
        error('mpc_constructors:WeightDimension', ...
            '%s must have %d elements.',name,dim);
    end
    W = diag(double(W(:)));
else
    if ~isequal(size(W),[dim dim])
        error('mpc_constructors:WeightDimension', ...
            '%s must be a %d-by-%d matrix.',name,dim,dim);
    end
    W = double(W);
end
if norm(W-W','fro') > 1e-10*max(1,norm(W,'fro'))
    error('mpc_constructors:WeightSymmetry','%s must be symmetric.',name);
end
W = (W + W')/2;
eigW = eig(W);
tol = 1e-10*max(1,max(abs(eigW)));
if any(eigW < -tol)
    error('mpc_constructors:WeightDefinite', ...
        '%s must be positive semidefinite.',name);
end
end

function s = localDefault(s,name,value)
if ~isfield(s,name) || isempty(s.(name))
    s.(name) = value;
end
end
