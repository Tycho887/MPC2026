function [estimate, filter] = problem06Filter(model, mode, y, u)
%PROBLEM06FILTER Design and run a static or dynamic Kalman filter.
%
%   [estimate, filter] = problem06Filter(model, mode, y, u)
%
% model must contain A, B, G, C, Q, R, x0 and P0. D and Ts are optional.
% mode is 'static' or 'dynamic'.
%
% The static design uses the Control System Toolbox KALMAN function. The
% dynamic design uses the equivalent time-varying Riccati recursion because
% KALMAN returns only the steady-state solution. Both filters use the same
% prediction/measurement-update bookkeeping, so their estimates are directly
% comparable.

mode = validatestring(mode,{'static','dynamic'},mfilename,'mode');

A = model.A;
B = model.B;
G = model.G;
C = model.C;
Q = model.Q;
R = model.R;
if isfield(model,'D') && ~isempty(model.D)
    D = model.D;
else
    D = zeros(size(C,1),size(B,2));
end
if isfield(model,'Ts') && ~isempty(model.Ts)
    Ts = model.Ts;
else
    Ts = 1;
end

nx = size(A,1);
ny = size(C,1);
ng = size(G,2);

if strcmp(mode,'static')
    sys = ss(A,[B G],C,[D zeros(ny,ng)],Ts);
    [kest,L,P,Mx] = kalman(sys,Q,R);
    if isempty(Mx)
        error('kalman did not return a discrete measurement-update gain.');
    end
    filter = struct('mode',mode,'L',L,'P',P,'Mx',Mx,'kest',kest);
else
    filter = struct('mode',mode,'L',[],'P',[],'Mx',[],'kest',[]);
end
filter.P0 = model.P0;
filter.x0 = model.x0;
filter.Ts = Ts;

if nargin < 3 || isempty(y) || nargin < 4 || isempty(u)
    estimate = struct();
    return;
end

N = size(y,2);
if size(u,2) ~= N
    error('y and u must have the same number of samples.');
end
if size(y,1) ~= ny || size(u,1) ~= size(B,2)
    error('y or u has incompatible dimensions.');
end

xPred = zeros(nx,N);
xFilt = zeros(nx,N);
innovation = zeros(ny,N);
PHist = zeros(nx,nx,N);
KHist = zeros(nx,ny,N);

xPred(:,1) = model.x0;
if strcmp(mode,'static')
    P = filter.P;
else
    P = model.P0;
end

for k = 1:N
    e = y(:,k) - C*xPred(:,k);

    if strcmp(mode,'static')
        K = filter.Mx;
        xFilt(:,k) = xPred(:,k) + K*e;
    else
        Re = C*P*C' + R;
        K = P*C'/Re;
        xFilt(:,k) = xPred(:,k) + K*e;
        Pf = P - K*Re*K';
        Pf = (Pf + Pf')/2;
    end

    innovation(:,k) = e;
    PHist(:,:,k) = P;
    KHist(:,:,k) = K;

    if k < N
        xPred(:,k+1) = A*xFilt(:,k) + B*u(:,k);
    end

    if strcmp(mode,'dynamic')
        P = A*Pf*A' + G*Q*G';
        P = (P + P')/2;
    end
end

filter.PFinal = P;
filter.finalGain = KHist(:,:,end);

estimate = struct();
estimate.mode = mode;
estimate.xPred = xPred;
estimate.xFilt = xFilt;
estimate.zPred = C*xPred;
estimate.zFilt = C*xFilt;
estimate.innovation = innovation;
estimate.P = PHist;
estimate.K = KHist;
end
