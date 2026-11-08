function model = problem06Model(config)
%PROBLEM06MODEL Discrete linear four-tank model for Problems 6.2 and 6.3.
% Returns the 4-state mass-deviation model and the 6-state disturbance-
% augmented model used by the offset-free filter.
%
%   model.base      A,B,G,C,D,Q,R,x0,P0,Ts  (4 states)
%   model.augmented A,B,G,C,D,Q,R,x0,P0,Ts  (6 states: masses + F3,F4)
%
% All states and inputs are deviations from the Problem 5 operating point.
% The process noise enters through the physical disturbance matrix E.

if nargin < 1 || isempty(config)
    config = struct();
end
config = localDefault(config,'Ts',5);
config = localDefault(config,'u0',[300;300]);
config = localDefault(config,'d0',[20;20]);
config = localDefault(config,'sigmaF',5);
config = localDefault(config,'sigmaV',0.5);
config = localDefault(config,'sigmaB',2);
config = localDefault(config,'P0',diag([100 100 50 50].^2));
config = localDefault(config,'P0b',diag([5 5].^2));

p = fourTankParameters();
op = fourTankSteadyState(config.u0,config.d0,p);
continuous = linearizeFourTank(op,p);
discrete = discretizeLinearModel(continuous,config.Ts,'zoh');

A = discrete.A;
B = discrete.B;
E = discrete.E;
C = discrete.C;
D = discrete.D;
Q = config.sigmaF^2 * eye(2);
R = config.sigmaV^2 * eye(2);

base = struct('A',A,'B',B,'G',E,'C',C,'D',D,'Q',Q,'R',R, ...
    'x0',zeros(4,1),'P0',config.P0,'Ts',config.Ts);

% The disturbance bias d_k enters the mass states through E and evolves as
% d_{k+1} = d_k + w_b,k. The white part of the disturbance enters through E
% in the same way as in the 4-state model.
Aa = [A, E; zeros(2,4), eye(2)];
Ba = [B; zeros(2,2)];
Ga = [E, zeros(4,2); zeros(2,2), eye(2)];
Ca = [C, zeros(2,2)];
Da = zeros(2,2);
Qa = blkdiag(Q,config.sigmaB^2 * eye(2));
P0a = blkdiag(config.P0,config.P0b);

augmented = struct('A',Aa,'B',Ba,'G',Ga,'C',Ca,'D',Da,'Q',Qa,'R',R, ...
    'x0',zeros(6,1),'P0',P0a,'Ts',config.Ts);

model = struct();
model.p = p;
model.op = op;
model.continuous = continuous;
model.discrete = discrete;
model.config = config;
model.base = base;
model.augmented = augmented;
model.stateNames = {'m1','m2','m3','m4'};
model.augmentedStateNames = {'m1','m2','m3','m4','F3','F4'};

assert(rank(obsv(A,C)) == size(A,1), ...
    'The 4-state mass model is not observable from the lower-tank levels.');
assert(rank(obsv(Aa,Ca)) == size(Aa,1), ...
    'The disturbance-augmented model is not observable from the lower-tank levels.');
end

function s = localDefault(s,name,value)
if ~isfield(s,name) || isempty(s.(name))
    s.(name) = value;
end
end
