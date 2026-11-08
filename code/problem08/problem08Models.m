function models = problem08Models(config)
%PROBLEM08MODELS Linear models used by the Problem 8 unconstrained MPC.
%
%   models = problem08Models(config)
%
% models.linear      Problem 5 mass model returned by problem06Model. The
%                    4-state design matrices A, B, G, C, D and Ts are also
%                    exposed at the top level.
% models.identified  6-state ZOH realization of the Problem 4 identified
%                    transfer-function matrix, or [] if unavailable. The
%                    original continuous model is kept in transferFunction,
%                    not G, so the struct is a valid MPC design model.
%
% Optional config fields: Ts (5 s), identificationFile, and
% identificationBackend ('least_squares').

if nargin < 1 || isempty(config)
    config = struct();
end
config = localDefault(config,'Ts',5);
config = localDefault(config,'identificationFile',fullfile( ...
    fileparts(fileparts(mfilename('fullpath'))), ...
    'problem04','results','identification.mat'));
config = localDefault(config,'identificationBackend','least_squares');
backend = validatestring(config.identificationBackend, ...
    {'least_squares','toolbox','auto'},mfilename,'identificationBackend');

models = struct();
models.linear = problem06Model(config);
base = models.linear.base;
models.linear.A = base.A;
models.linear.B = base.B;
models.linear.G = base.G;
models.linear.C = base.C;
models.linear.D = base.D;
models.linear.Ts = base.Ts;

models.identified = localIdentified(config.identificationFile,backend, ...
    config.Ts);
end

function identified = localIdentified(file,backend,Ts)
identified = [];
try
    if isfile(file)
        [G,op] = localLoadIdentification(file);
    else
        [G,op] = localComputeIdentification(backend);
    end
    identified = localDiscretize(G,op,Ts);
catch err
    warning('problem08Models:identified', ...
        'Identified model unavailable: %s',err.message);
end
end

function [G,op] = localLoadIdentification(file)
saved = load(file);
if isfield(saved,'result')
    saved = saved.result;
end
if ~isfield(saved,'G')
    error('No identified transfer function in %s.',file);
end
G = saved.G;
op = struct();
if isfield(saved,'op')
    op = saved.op;
end
end

function [G,op] = localComputeIdentification(backend)
problem04Directory = fullfile( ...
    fileparts(fileparts(mfilename('fullpath'))),'problem04');
onPath = any(strcmp(strsplit(path,pathsep),problem04Directory));
if ~onPath
    addpath(problem04Directory);
end
restorePath = onCleanup(@() localRestorePath(problem04Directory,onPath));
[result, figures] = problem04Identification(backend);
close(figures(isgraphics(figures)));
G = result.G;
op = result.op;
end

function localRestorePath(directory,onPath)
if ~onPath && any(strcmp(strsplit(path,pathsep),directory))
    rmpath(directory);
end
end

function identified = localDiscretize(G,op,Ts)
discrete = c2d(G,Ts,'zoh');
sys = ss(discrete);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;
assert(norm(D,1) <= 1e-10, ...
    'The identified discrete model must have zero feedthrough.');
D = zeros(size(D));
identified = struct('A',A,'B',B,'C',C,'D',D,'Ts',Ts,'op',op, ...
    'transferFunction',G);
end

function s = localDefault(s,name,value)
if ~isfield(s,name) || isempty(s.(name))
    s.(name) = value;
end
end
