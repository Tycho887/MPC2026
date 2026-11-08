function estimator = mpc_estimators(kind, model, mode)
%MPC_ESTIMATORS Build a stationary Kalman estimator struct for the MPC.
%
%   estimator = mpc_estimators(kind, model, mode)
%
% kind is 'unbiased' (4-state mass model) or 'biased' (6-state disturbance-
% augmented model). mode is 'static' (default). model is a problem06Model
% struct; when omitted it is built with problem06Model.
%
% The MPC uses the steady-state gains, so mode='dynamic' is rejected here.
% Dynamic filters remain available in problem06Filter for the Problem 6 study.

if nargin < 1 || isempty(kind)
    kind = 'unbiased';
end
kind = validatestring(kind,{'unbiased','biased'},mfilename,'kind');
if nargin < 2 || isempty(model)
    model = problem06Model();
end
if nargin < 3 || isempty(mode)
    mode = 'static';
end
mode = validatestring(mode,{'static','dynamic'},mfilename,'mode');
if strcmp(mode,'dynamic')
    error('mpc_estimators:dynamicMode', ...
        ['mode=''dynamic'' is not supported: the MPC uses stationary filter ' ...
         'gains. Use problem06Filter for the dynamic Problem 6 filters.']);
end

assert(isstruct(model) && isfield(model,'base'), ...
    'mpc_estimators:model','model must be a problem06Model struct with a base field.');
base = model.base;
nx = size(base.A,1);
assert(size(base.A,2) == nx && size(base.C,2) == nx, ...
    'mpc_estimators:model','base A and C must describe the same state dimension.');

config = struct();
if isfield(model,'config') && ~isempty(model.config)
    config = model.config;
end

switch kind
    case 'unbiased'
        [~, filter] = problem06Filter(base,'static');
        nw = 0;
        Kfx = filter.Mx;
        Kfw = [];
        isPersistent = false;
        x0 = base.x0;
        w0 = [];
        design = base;
    case 'biased'
        assert(isfield(model,'augmented') && ~isempty(model.augmented), ...
            'mpc_estimators:model','kind=''biased'' needs model.augmented.');
        aug = model.augmented;
        [~, filter] = problem06Filter(aug,'static');
        nw = size(aug.A,1) - nx;
        assert(nw > 0 && size(aug.A,1) == size(aug.C,2), ...
            'mpc_estimators:model','augmented A and C dimensions are inconsistent.');
        Kfx = filter.Mx(1:nx,:);
        Kfw = filter.Mx(nx+1:nx+nw,:);
        isPersistent = true;
        x0 = aug.x0(1:nx);
        w0 = aug.x0(nx+1:nx+nw);
        design = aug;
end

ny = size(design.C,1);
assert(isequal(size(Kfx),[nx ny]), ...
    'mpc_estimators:gain','Kfx must be %d-by-%d, got %d-by-%d.', ...
    nx,ny,size(Kfx,1),size(Kfx,2));
if nw > 0
    assert(isequal(size(Kfw),[nw ny]), ...
        'mpc_estimators:gain','Kfw must be %d-by-%d, got %d-by-%d.', ...
        nw,ny,size(Kfw,1),size(Kfw,2));
end

sigmaF = localValue(config,'sigmaF',sqrt(base.Q(1,1)));
sigmaV = localValue(config,'sigmaV',sqrt(base.R(1,1)));
if isfield(model,'augmented') && ~isempty(model.augmented) && nx < size(model.augmented.A,1)
    nq = size(base.Q,1);
    sigmaB = localValue(config,'sigmaB',sqrt(model.augmented.Q(nq+1,nq+1)));
else
    sigmaB = localValue(config,'sigmaB',NaN);
end

estimator = struct();
estimator.kind = kind;
estimator.mode = mode;
estimator.nx = nx;
estimator.nw = nw;
estimator.Kfx = Kfx;
estimator.Kfw = Kfw;
estimator.persistent = isPersistent;
estimator.x0 = x0;
estimator.w0 = w0;
estimator.model = design;
estimator.filter = filter;
estimator.sigmaF = sigmaF;
estimator.sigmaV = sigmaV;
estimator.sigmaB = sigmaB;
end

function value = localValue(s,name,fallback)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value = s.(name);
else
    value = fallback;
end
end
