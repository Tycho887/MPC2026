function validateMigration(includeEntryScripts)
%VALIDATEMIGRATION Numerical checks; optional full report artifact smoke runs.
if nargin < 1, includeEntryScripts = false; end
root = fileparts(fileparts(mfilename('fullpath')));
oldDirectory = pwd; oldPath = path;
visibility = get(groot,'DefaultFigureVisible');
existingFigures = findall(groot,'Type','figure');
restore = onCleanup(@() restoreSession(oldDirectory,oldPath,visibility,existingFigures)); %#ok<NASGU>
addpath(fullfile(root,'common'),fullfile(root,'problem02'), ...
    fullfile(root,'problem03'),fullfile(root,'problem04'));
set(groot,'DefaultFigureVisible','off');
p = fourTankParameters();
operatingInputs = {[300;300],[100;300],[300;300],[300;300]};
operatingDisturbances = {[20;20],[100;100],[0;0],[100;100]};
for j = 1:4
    op = fourTankSteadyState(operatingInputs{j},operatingDisturbances{j},p);
    assert(norm(fourTankDynamics(0,op.x,op.u,op.d,p),inf)<1e-10);
    model = linearizeFourTank(op,p);
    numericalA = zeros(4); numericalB = zeros(4,2); numericalE = zeros(4,2);
    for i = 1:4
        delta = zeros(4,1); delta(i) = 1e-4*op.x(i);
        numericalA(:,i) = (fourTankDynamics(0,op.x+delta,op.u,op.d,p) ...
            - fourTankDynamics(0,op.x-delta,op.u,op.d,p))/(2*delta(i));
    end
    for i = 1:2
        delta = zeros(2,1); delta(i) = 0.01;
        numericalB(:,i) = (fourTankDynamics(0,op.x,op.u+delta,op.d,p) ...
            - fourTankDynamics(0,op.x,op.u-delta,op.d,p))/(2*delta(i));
        numericalE(:,i) = (fourTankDynamics(0,op.x,op.u,op.d+delta,p) ...
            - fourTankDynamics(0,op.x,op.u,op.d-delta,p))/(2*delta(i));
    end
    assert(norm(model.A-numericalA,inf)<1e-8);
    assert(norm(model.B-numericalB,inf)<1e-9);
    assert(norm(model.E-numericalE,inf)<1e-9);
    assert(norm(fourTankLevels(op.x,p)-op.h,inf)<1e-12);
end
% Nonuniform tank areas and density exercise the mass/level scaling.
pUnequal = p; pUnequal.A = p.A.*[1;2;3;4]; pUnequal.rho = 1.2;
op = fourTankSteadyState([300;300],[20;20],pUnequal);
model = linearizeFourTank(op,pUnequal);
[Gm,Gh] = fourTankTransferFunctions(model);
assert(norm(dcgain(Gh)-diag(1./(pUnequal.rho*pUnequal.A(1:2)))*dcgain(Gm),inf)<1e-8);
controller = tuneFourTankPI(Gh);
assert(norm(dcgain(Gh)*controller.W-eye(2),inf)<1e-10);

[r2,~] = problem02Experiments();
assert(isequal(r2.deterministic.y,r2.deterministic.z));
assert(isequal(size(r2.sde.x),[4 numel(r2.sde.t)]));
assert(all(r2.sde.x(:)>=0) && all(isfinite(r2.sde.x(:))));
[r2Repeat,~] = problem02Experiments();
assert(isequal(r2.sde.x,r2Repeat.sde.x));
assert(isequal(r2.piecewise_stochastic.y,r2Repeat.piecewise_stochastic.y));
assert(any(r2.piecewise_stochastic.d(:)~=20));
% Both ODE backends keep an equilibrium at rest.
for solver = {'ode15s','ode45'}
    r = simulateFourTank([0;10;20],op.x,repmat(op.u,1,3), ...
        repmat(op.d,1,3),pUnequal,struct('solver',solver{1}));
    assert(max(abs(r.x-op.x),[],'all')<1e-6);
end
[comparison,~] = problem03Comparison();
for j = 1:numel(comparison.experiments)
    e = comparison.experiments{j};
    assert(all(isfinite(e.xNonlinear(:))) && all(isfinite(e.xLinear(:))));
    assert(all(e.uNonlinear(:)>=0 & e.uNonlinear(:)<=1000));
    assert(all(e.uLinear(:)>=0 & e.uLinear(:)<=1000));
    assert(norm(e.xLinear-(e.xLinearDeviation+comparison.op.x),inf)<1e-10);
    assert(norm(e.uNonlinear-e.uLinear,inf)>1e-3);
end
[level,~] = problem03LevelControl();
assert(all(isfinite(level.x(:))) && all(level.u(:)>=0));
assert(max(abs(level.z(:,end)-level.reference(:,end)))<0.1);

[identifiedFallback,~] = problem04Identification('least_squares');
assert(strcmp(identifiedFallback.backend,'least_squares'));
assert(all(isfinite(dcgain(identifiedFallback.G)),'all'));
checkMarkov(identifiedFallback.G);
if exist('tfest','file') && license('test','Identification_Toolbox')
    [identifiedToolbox,~] = problem04Identification('toolbox');
    assert(strcmp(identifiedToolbox.backend,'toolbox'));
    checkMarkov(identifiedToolbox.G);
else
    fprintf('System Identification Toolbox unavailable: toolbox route skipped.\n');
end

source = fullfile(root,'..','..','MPC','four-tank-model');
if isfolder(source)
    reference = legacyResult(source,'problem2_4_3.m', ...
        "struct('x',x,'y',y,'z',z)");
    assert(max(abs(reference.x-r2.sde.x),[],'all')<1e-8);
    assert(max(abs(reference.y-r2.sde.y),[],'all')<1e-8);
    reference = legacyResult(source,'problem3_2.m', ...
        "struct('y',y_history','u',u_history')");
    assert(max(abs(reference.y-level.y),[],'all')<1e-4);
    assert(max(abs(reference.u-level.u),[],'all')<1e-3);
    fprintf('Legacy SDE and nonlinear level-control trajectories agree.\n');
end
if includeEntryScripts
    for local = [false true]
        for number = 2:4
            directory = fullfile(root,sprintf('problem%02d',number));
            figuresBefore = findall(groot,'Type','figure');
            if local
                cd(directory); eval('simulate');
            else
                cd(root); run(fullfile(directory,'simulate.m'));
            end
            current = findall(groot,'Type','figure');
            created = current(~ismember(current,figuresBefore));
            for i = 1:numel(created)
                axesHandles = findall(created(i),'Type','axes');
                for j = 1:numel(axesHandles)
                    assert(norm(axesHandles(j).Color-[1 1 1])<1e-10);
                    assert(norm(axesHandles(j).XColor)<1);
                end
            end
            mats = dir(fullfile(directory,'results','*.mat'));
            pngs = dir(fullfile(directory,'results','*.png'));
            pdfs = dir(fullfile(directory,'results','*.pdf'));
            expectedMats = [1 3 4]; expectedFigures = [3 5 6];
            assert(numel(mats)==expectedMats(number-1));
            assert(numel(pngs)==expectedFigures(number-1) && numel(pdfs)==numel(pngs));
            for i = 1:numel(mats)
                saved = load(fullfile(mats(i).folder,mats(i).name));
                assert(isfield(saved,'result'));
            end
            closeCreatedFigures(existingFigures);
        end
    end
    fprintf('All entry scripts and MAT/PNG/PDF exports passed from root and local directories.\n');
end
fprintf('Migration numerical checks passed.\n');
end

function checkMarkov(G)
r = markovParameters(G,4,60);
[stepResponse,~] = step(r.system,r.t);
% step returns time-by-output-by-input; coefficients use output-by-input-by-time.
S = permute(stepResponse,[2 3 1]);
expected = cat(3,S(:,:,1),diff(S,1,3));
assert(max(abs(expected-r.H),[],'all')<1e-8);
end

function result = legacyResult(source,filename,expression)
% Temporarily wrap a source script to capture its results without editing it.
directory = tempname; mkdir(directory);
cleanup = onCleanup(@() removeReference(directory)); %#ok<NASGU>
text = fileread(fullfile(source,filename));
text = strrep(text,'clear; clc; close all;','');
firstFunction = regexp(text,'(?m)^function ','once');
if isempty(firstFunction), firstFunction = numel(text)+1; end
body = text(1:firstFunction-1); helpers = text(firstFunction:end);
wrapper = sprintf('function result = migrationLegacyReference()\n%s\nresult = %s;\nend\n%s\n',body,expression,helpers);
fid = fopen(fullfile(directory,'migrationLegacyReference.m'),'w');
assert(fid>=0); fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',wrapper); clear fileCleanup;
addpath(directory); result = migrationLegacyReference();
end

function removeReference(directory)
rmpath(directory); clear migrationLegacyReference;
rmdir(directory,'s');
end

function restoreSession(directory,oldPath,visibility,figures)
cd(directory); path(oldPath); set(groot,'DefaultFigureVisible',visibility);
closeCreatedFigures(figures);
end

function closeCreatedFigures(existing)
current = findall(groot,'Type','figure');
created = current(~ismember(current,existing));
close(created);
end
