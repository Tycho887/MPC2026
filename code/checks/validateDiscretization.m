function validateDiscretization(includeEntryScript)
%VALIDATEDISCRETIZATION Numerical checks for Problem 5 and shared Markov code.
% Pass true to also run Problem 5 from root and locally, checking exports.
if nargin < 1, includeEntryScript = false; end
root = fileparts(fileparts(mfilename('fullpath')));
oldDirectory = pwd; oldPath = path;
visibility = get(groot,'DefaultFigureVisible');
existingFigures = findall(groot,'Type','figure');
cleanup = onCleanup(@() restoreSession(oldDirectory,oldPath,visibility,existingFigures));
addpath(fullfile(root,'common'));
set(groot,'DefaultFigureVisible','off');
p = fourTankParameters();
op = fourTankSteadyState([300;300],[20;20],p);
continuous = linearizeFourTank(op,p);
continuous.Cz = continuous.C;

for Ts = [0.5 4]
    zoh = discretizeLinearModel(continuous,Ts,'zoh');
    exponential = discretizeLinearModel(continuous,Ts,'expm');
    for field = {'A','B','E'}
        assert(max(abs(zoh.(field{1})-exponential.(field{1})),[],'all')<1e-10);
    end
    assert(isequal(zoh.op,op) && isequal(exponential.Cz,continuous.C));
    Hzoh = markovParameters(ss(zoh.A,zoh.B,zoh.C,zoh.D,Ts),Ts,60);
    Hexpm = markovParameters(ss(exponential.A,exponential.B, ...
        exponential.C,exponential.D,Ts),Ts,60);
    assert(isequal(size(Hzoh.H),[2 2 60]));
    assert(max(abs(Hzoh.H-Hexpm.H),[],'all')<1e-10);
    assert(isequal(Hzoh.H(:,:,1),zoh.D));
    assert(norm(Hzoh.H(:,:,2)-zoh.C*zoh.B,inf)<1e-12);
    assert(norm(Hzoh.H(:,:,3)-zoh.C*zoh.A*zoh.B,inf)<1e-12);
    checkStep(Hzoh);
    checkStep(Hexpm);
end

% A nilpotent integrator has singular A and closed-form ZOH matrices.
Ts = 4;
integrator = struct('A',[0 1;0 0],'B',[0;1],'E',[1;0], ...
    'C',eye(2),'D',zeros(2,1));
for method = {'zoh','expm'}
    discrete = discretizeLinearModel(integrator,Ts,method{1});
    assert(norm(discrete.A-[1 Ts;0 1],inf)<1e-12);
    assert(norm(discrete.B-[Ts^2/2;Ts],inf)<1e-12);
    assert(norm(discrete.E-[Ts;0],inf)<1e-12);
end

% Discrete input keeps an unobservable state and includes nonzero H_0.
system = ss(blkdiag(discrete.A,0.5),[discrete.B;0], ...
    [discrete.C zeros(2,1)],[0.25;0.5],Ts);
markov = markovParameters(system,Ts,8);
assert(isequal(markov.system.A,system.A));
assert(isequal(markov.H(:,:,1),system.D));
checkStep(markov);
single = markovParameters(system,Ts,1);
assert(isequal(single.H,system.D) && isequal(single.t,0));
try
    markovParameters(system,2,8);
    error('validateDiscretization:MissingError','A mismatched sample time must fail.');
catch exception
    assert(strcmp(exception.identifier,'markovParameters:SamplingTimeMismatch'));
end

% Problem 4's continuous transfer-function path retains its previous behavior.
[~,G] = fourTankTransferFunctions(continuous);
markov = markovParameters(G,Ts,60);
previous = minreal(ss(c2d(G,Ts,'zoh')),1e-8);
assert(isequal(markov.system.A,previous.A));
assert(isequal(markov.system.B,previous.B));
checkStep(markov);

if includeEntryScript
    directory = fullfile(root,'problem05');
    for local = [false true]
        if local
            cd(directory); problem5_4;
        else
            cd(fileparts(root)); run(fullfile(directory,'problem5_4.m'));
        end
        assert(problem05Results.differences.H<1e-10);
        saved = load(fullfile(directory,'results','discretization.mat'));
        assert(isequal(saved.result.zoh.markov.H,problem05Results.zoh.markov.H));
        for extension = {'png','pdf'}
            artifact = dir(fullfile(directory,'results',['discretization_01.' extension{1}]));
            assert(isscalar(artifact) && artifact.bytes>0);
        end
        assert(numel(findall(problem05Figures,'Type','axes'))==4);
        close(problem05Figures);
    end
    fprintf('Problem 5 entry script and MAT/PNG/PDF exports passed from root and locally.\n');
end
fprintf('Discretization and Markov numerical checks passed.\n');
end

function checkStep(markov)
[response,~] = step(markov.system,markov.t);
expected = permute(response,[2 3 1]);
assert(max(abs(cumsum(markov.H,3)-expected),[],'all')<1e-10);
end

function restoreSession(directory,oldPath,visibility,existingFigures)
current = findall(groot,'Type','figure');
close(current(~ismember(current,existingFigures)));
cd(directory); path(oldPath);
set(groot,'DefaultFigureVisible',visibility);
end
