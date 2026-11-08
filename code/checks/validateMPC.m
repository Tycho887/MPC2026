function validateMPC(includeEntryScript)
%VALIDATEMPC Numerical checks for the Problem 8 unconstrained MPC.
% Pass true to also run Problem 8 from root and locally, checking exports.
if nargin < 1, includeEntryScript = false; end
root = fileparts(fileparts(mfilename('fullpath')));
oldDirectory = pwd; oldPath = path;
visibility = get(groot,'DefaultFigureVisible');
existingFigures = findall(groot,'Type','figure');
cleanup = onCleanup(@() restoreSession(oldDirectory,oldPath,visibility, ...
    existingFigures));
addpath(fullfile(root,'common'),fullfile(root,'problem06'), ...
    fullfile(root,'problem07'),fullfile(root,'problem08'));
set(groot,'DefaultFigureVisible','off');

%% Random design model: prediction matrices and cost
fprintf('Checking prediction and cost matrices on a random 3-state model...\n');
model = localRandomModel();
N = 6;
config = struct('N',N,'Wz',diag([1.2 0.8]),'Wu',0.3*eye(2), ...
    'Wdu',0.5*eye(2));
rng(802,'twister');
estimator = struct('Kfx',0.2*randn(3,2),'Kfw',0.1*randn(2,2), ...
    'nw',2,'persistent',true);
mpc = mpc_constructors(model,config,estimator);
nx = size(model.A,1); nu = size(model.B,2);
nz = size(model.C,1); nw = size(model.G,2);

assert(isequal(size(mpc.Gamma),[nz*N nu*N]));
assert(isequal(size(mpc.Phi_x),[nz*N nx]));
assert(isequal(size(mpc.Phi_w),[nz*N nw]));
assert(mpc.nw == nw && mpc.persistent && ...
    strcmp(mpc.disturbanceMode,'persistent'));

rng(803,'twister');
x = randn(nx,1); w = randn(nw,1); U = randn(nu*N,1);
combined = mpc.Phi_x*x + mpc.Phi_w*w + mpc.Gamma*U;
assert(norm(localRollout(model,x,w,U,N)-combined,inf) < 1e-10);
assert(norm(localRollout(model,zeros(nx,1),zeros(nw,1),U,N) - ...
    mpc.Gamma*U,inf) < 1e-10);
assert(norm(localRollout(model,x,zeros(nw,1),zeros(nu*N,1),N) - ...
    mpc.Phi_x*x,inf) < 1e-10);
assert(norm(localRollout(model,zeros(nx,1),w,zeros(nu*N,1),N) - ...
    mpc.Phi_w*w,inf) < 1e-10);

Wz_bar = kron(eye(N),mpc.Wz);
Wu_bar = kron(eye(N),mpc.Wu);
Wdu_bar = kron(eye(N),mpc.Wdu);
assert(norm(mpc.Wz_bar-Wz_bar,'fro') == 0);
assert(norm(mpc.Wu_bar-Wu_bar,'fro') == 0);
assert(norm(mpc.Wdu_bar-Wdu_bar,'fro') == 0);

Utest = reshape(1:nu*N,nu,N);
uTestPrev = (1:nu)';
delta = mpc.Lambda*Utest(:) - mpc.I0*uTestPrev;
deltaExpected = [Utest(:,1)-uTestPrev; reshape(diff(Utest,1,2),nu*(N-1),1)];
assert(norm(delta-deltaExpected,inf) == 0);

Hexpected = mpc.Gamma'*(Wz_bar'*Wz_bar)*mpc.Gamma + Wu_bar'*Wu_bar + ...
    (Wdu_bar*mpc.Lambda)'*(Wdu_bar*mpc.Lambda);
assert(norm(mpc.H-Hexpected,'fro') < 1e-10*max(1,norm(Hexpected,'fro')));
assert(norm(mpc.H-mpc.H','fro') < 1e-12*max(1,norm(mpc.H,'fro')));
assert(all(eig(mpc.H) > 0));
Mz = -mpc.Gamma'*(Wz_bar'*Wz_bar);
Mu = -(Wu_bar'*Wu_bar);
Mdu = -(Wdu_bar*mpc.Lambda)'*(Wdu_bar*mpc.I0);
assert(norm(mpc.Mz-Mz,'fro') < 1e-10*max(1,norm(Mz,'fro')));
assert(norm(mpc.Mu-Mu,'fro') < 1e-10*max(1,norm(Mu,'fro')));
assert(norm(mpc.Mdu-Mdu,'fro') < 1e-10*max(1,norm(Mdu,'fro')));

%% Unconstrained solve
fprintf('Checking the unconstrained solve and quadprog agreement...\n');
rng(804,'twister');
xPrev = randn(nx,1); wPrev = randn(nw,1); uPrev = randn(nu,1);
y = randn(nz,1); R = randn(nz*N,1); U_bar = repmat(randn(nu,1),N,1);
[uk,xHat,wHat,info] = mpc_step(mpc,xPrev,wPrev,uPrev,y,R,U_bar);
assert(isequal(uk,info.U(1:nu)));
assert(isfield(info,'solverInfo') && info.solverInfo == 1);
assert(isequal(size(info.qp.A),[0 nu*N]));
assert(all(isinf(info.qp.l)) && all(info.qp.l < 0));
assert(all(isinf(info.qp.u)) && all(info.qp.u > 0));
assert(isequal(info.qp.H,mpc.H));
assert(norm(mpc.H*info.U+info.g,inf) < 1e-8);
xPred = model.A*xPrev + model.B*uPrev + model.G*wPrev;
assert(norm(info.innovation-(y-model.C*xPred),inf) < 1e-12);
b = mpc.Phi_x*xHat + mpc.Phi_w*wHat;
assert(norm(info.b-b,inf) < 1e-12);
assert(norm(info.c-(R-b),inf) < 1e-12);
assert(norm(info.Z(:)-(b+mpc.Gamma*info.U),inf) < 1e-9);
options = optimoptions('quadprog','Display','off');
Uquad = quadprog(mpc.H,info.g,[],[],[],[],[],[],[],options);
assert(norm(info.U-Uquad,inf) < 1e-6);

%% Horizon rollout
fprintf('Checking the horizon rollout...\n');
rng(805,'twister');
x = randn(nx,1); w = randn(nw,1); U = randn(nu*N,1);
prediction = mpc_predict(mpc,x,w,U);
assert(norm(prediction.Z(:)-(prediction.b+mpc.Gamma*U),inf) < 1e-10);
assert(norm(prediction.Z(:)-localRollout(model,x,w,U,N),inf) < 1e-10);
assert(isequal(prediction.X(:,1),x));
assert(norm(prediction.U-reshape(U,nu,N),inf) == 0);

%% Steady-state target
fprintf('Checking steady-state targets...\n');
rng(806,'twister');
r = randn(nz,1); w = randn(nw,1);
[xss,uss] = mpc_target(mpc,r,w);
assert(norm(model.C*xss-r,inf) < 1e-10);
assert(norm(xss-(model.A*xss+model.B*uss+model.G*w),inf) < 1e-10);

%% Estimators against problem06Filter
fprintf('Checking both estimators against problem06Filter...\n');
models = problem06Model();
estUnbiased = mpc_estimators('unbiased',models);
estBiased = mpc_estimators('biased',models);
assert(estUnbiased.nw == 0 && ~estUnbiased.persistent && ...
    isempty(estUnbiased.Kfw));
assert(estBiased.nw == 2 && estBiased.persistent);
mpcUnbiased = mpc_constructors(models.base,struct('N',N),estUnbiased);
mpcBiased = mpc_constructors(models.base,struct('N',N),estBiased);

rng(807,'twister');
samples = 8;
y = randn(2,samples);
u = 10*randn(2,samples);
[estimate,filterUnbiased] = problem06Filter(models.base,'static',y,u);
[estimateAugmented,filterAugmented] = ...
    problem06Filter(models.augmented,'static',y,u);
assert(norm(estUnbiased.Kfx-filterUnbiased.Mx,inf) < 1e-12);
assert(norm(estBiased.Kfx-filterAugmented.Mx(1:4,:),inf) < 1e-12);
assert(norm(estBiased.Kfw-filterAugmented.Mx(5:6,:),inf) < 1e-12);

xHatPrev = estimate.xFilt(:,1);
for k = 2:samples
    [xHat,wHat] = mpc_estimate(mpcUnbiased,xHatPrev,[],u(:,k-1),y(:,k));
    assert(isempty(wHat));
    assert(norm(xHat-estimate.xFilt(:,k),inf) < 1e-10);
    xHatPrev = xHat;
end
augPrev = estimateAugmented.xFilt(:,1);
for k = 2:samples
    [xHat,wHat] = mpc_estimate(mpcBiased,augPrev(1:4),augPrev(5:6), ...
        u(:,k-1),y(:,k));
    augPrev = [xHat; wHat];
    assert(norm(augPrev-estimateAugmented.xFilt(:,k),inf) < 1e-10);
end

%% Deterministic offset contrast
fprintf(['Checking the deterministic offset contrast ' ...
    '(sigmaF = sigmaV = 1e-3)...\n']);
[result,figures] = problem08Simulate(struct('sigmaF',1e-3,'sigmaV',1e-3));
biasedOffset = result.metrics.biased.finalOffsetMaxAbs;
unbiasedOffset = result.metrics.unbiased.finalOffsetMaxAbs;
assert(biasedOffset < 0.1,'validateMPC:BiasedOffset', ...
    'The biased MPC must remove the offset, got %.4f cm.',biasedOffset);
assert(unbiasedOffset > 1,'validateMPC:UnbiasedOffset', ...
    'The unbiased MPC must keep an offset, got %.4f cm.',unbiasedOffset);
assert(numel(figures) == 4 && all(isgraphics(figures)), ...
    'validateMPC:SimulationFigures', ...
    'problem08Simulate must return four valid figures.');
close(figures);
fprintf('  biased offset %.4f cm, unbiased offset %.3f cm\n', ...
    biasedOffset,unbiasedOffset);

%% Entry script and artifacts
if includeEntryScript
    directory = fullfile(root,'problem08');
    for local = [false true]
        clear('problem08Result','problem08Figures');
        if local
            cd(directory); run_mpc;
        else
            cd(fileparts(root)); run(fullfile(directory,'run_mpc.m'));
        end
        assert(exist('problem08Result','var') == 1, ...
            'validateMPC:EntryResult','run_mpc must define problem08Result.');
        assert(numel(problem08Figures) == 4 && ...
            all(isgraphics(problem08Figures)), ...
            'validateMPC:EntryFigures', ...
            'run_mpc must leave four valid figures.');
        for k = 1:numel(problem08Figures)
            assert(~isempty(findall(problem08Figures(k),'Type','axes')), ...
                'validateMPC:EntryAxes', ...
                'Figure %d of run_mpc has no axes.',k);
        end

        savedFile = fullfile(directory,'results','unconstrained_mpc.mat');
        assert(isfile(savedFile),'validateMPC:EntryMat', ...
            'Missing %s.',savedFile);
        fileInfo = dir(savedFile);
        assert(fileInfo.bytes > 0);
        saved = load(savedFile);
        assert(isfield(saved,'result'));
        assert(abs(saved.result.metrics.biased.finalOffsetMaxAbs - ...
            problem08Result.metrics.biased.finalOffsetMaxAbs) < 1e-10);

        for k = 1:4
            for extension = {'png','pdf'}
                artifact = dir(fullfile(directory,'results', ...
                    sprintf('unconstrained_mpc_%02d.%s',k,extension{1})));
                assert(isscalar(artifact) && artifact.bytes > 0, ...
                    'validateMPC:EntryArtifact', ...
                    'Missing or empty unconstrained_mpc_%02d.%s.', ...
                    k,extension{1});
            end
        end
        close(problem08Figures);
    end
    fprintf(['Problem 8 entry script and MAT/PNG/PDF exports passed ' ...
        'from root and locally.\n']);
end

fprintf('MPC numerical checks passed.\n');
end

function model = localRandomModel()
rng(801,'twister');
nx = 3; nu = 2; nz = 2; nw = 2;
A = randn(nx);
A = 0.85*A/max(abs(eig(A)));
B = randn(nx,nu);
C = randn(nz,nx);
G = randn(nx,nw);
model = struct('A',A,'B',B,'C',C,'G',G,'D',zeros(nz,nu),'Ts',1);
end

function Z = localRollout(model,x,w,U,N)
nu = size(model.B,2);
Umat = reshape(U,nu,N);
nx = size(model.A,1);
X = zeros(nx,N+1);
X(:,1) = x;
Z = zeros(size(model.C,1),N);
for j = 1:N
    X(:,j+1) = model.A*X(:,j) + model.B*Umat(:,j) + model.G*w;
    Z(:,j) = model.C*X(:,j+1);
end
Z = Z(:);
end

function restoreSession(directory,oldPath,visibility,existingFigures)
current = findall(groot,'Type','figure');
close(current(~ismember(current,existingFigures)));
cd(directory); path(oldPath);
set(groot,'DefaultFigureVisible',visibility);
end
