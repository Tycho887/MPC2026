function [result, figures] = problem08Simulate(config)
%PROBLEM08SIMULATE Unconstrained MPC on the linear four-tank model.
%
%   [result, figures] = problem08Simulate(config)
%
% Compares two offset-free MPC designs on the same linear plant and noise
% realization:
%   unbiased  - 4-state mass model with no disturbance state
%   biased    - 6-state model with F3/F4 bias states
%
% The experiment applies a level reference step and two disturbance steps.
% It reports tracking metrics, plots the closed-loop responses, the bias
% estimates, the physical pump flows and one horizon prediction. A design-
% only MPC is built for the Problem 4 identified model when it is available.
% No files are written here; run_mpc.m saves the result.

if nargin < 1 || isempty(config)
    config = struct();
end
config = localDefault(config,'Ts',5);
config = localDefault(config,'u0',[300;300]);
config = localDefault(config,'d0',[20;20]);
config = localDefault(config,'N',20);
config = localDefault(config,'Wz',10*eye(2));
config = localDefault(config,'Wu',0.01*eye(2));
config = localDefault(config,'Wdu',1*eye(2));
config = localDefault(config,'tFinal',1000);
config = localDefault(config,'seed',8001);
config = localDefault(config,'refTime',200);
config = localDefault(config,'ref',[1;0.5]);
config = localDefault(config,'stepTimes',[300 500]);
config = localDefault(config,'stepSizes',[20;15]);
config = localDefault(config,'settleTime',100);
config = localDefault(config,'predictionTime',700);

%% Models, estimators and controllers
models = problem08Models(config);
base = models.linear.base;
op = models.linear.op;

estUnbiased = mpc_estimators('unbiased',models.linear);
estBiased = mpc_estimators('biased',models.linear);
mpcUnbiased = mpc_constructors(models.linear,config,estUnbiased);
mpcBiased = mpc_constructors(models.linear,config,estBiased);

mpcIdentified = [];
identifiedAvailable = ~isempty(models.identified);
if identifiedAvailable
    mpcIdentified = mpc_constructors(models.identified,config,[]);
end

%% Truth signals and one shared noise realization
Ts = base.Ts;
t = (0:Ts:config.tFinal)';
Nsim = numel(t);

sigmaF = models.linear.config.sigmaF;
sigmaV = models.linear.config.sigmaV;
rng(config.seed,'twister');
w = sigmaF*randn(2,Nsim);
v = sigmaV*randn(2,Nsim);

dTrue = zeros(2,Nsim);
stepTimes = config.stepTimes(:)';
stepSizes = config.stepSizes(:)';
for j = 1:numel(stepTimes)
    dTrue(j,t >= stepTimes(j)) = stepSizes(j);
end

r = zeros(2,Nsim);
r(:,t >= config.refTime) = repmat(config.ref(:),1,sum(t >= config.refTime));

%% Closed-loop simulation of both controllers
predictionTime = min(max(config.predictionTime,0),config.tFinal);
[~,predictionTarget] = min(abs(t - predictionTime));
predictionIndex = min(predictionTarget,max(1,Nsim - mpcBiased.N));

traceUnbiased = localSimulate(mpcUnbiased,base,dTrue,w,v,r,t,predictionIndex);
traceBiased = localSimulate(mpcBiased,base,dTrue,w,v,r,t,predictionIndex);

%% Metrics
metricsUnbiased = localMetrics(traceUnbiased,r,dTrue,config);
metricsBiased = localMetrics(traceBiased,r,dTrue,config);

%% Figures
hMeasScale = max(1,round(Nsim/60));

figTracking = experimentFigure('Name','Problem 8: reference tracking', ...
    'Position',[100 100 1200 800]);
sgtitle('Problem 8: reference tracking, unbiased (left) vs biased (right)');
traces = {traceUnbiased,traceBiased};
traceNames = {'unbiased','biased'};
axTrack = gobjects(2,2);
for col = 1:2
    trace = traces{col};
    for channel = 1:2
        ax = subplot(2,2,(channel-1)*2 + col);
        axTrack(channel,col) = ax;
        plot(ax,t,op.h(channel) + trace.z(channel,:),'k-', ...
            'LineWidth',1.3,'DisplayName','true z');
        hold(ax,'on');
        stairs(ax,t,op.h(channel) + r(channel,:),'b--', ...
            'LineWidth',1.2,'DisplayName','reference');
        idx = 1:hMeasScale:Nsim;
        plot(ax,t(idx),op.h(channel) + trace.y(channel,idx),'.', ...
            'Color',[0.45 0.45 0.45],'MarkerSize',7, ...
            'DisplayName','measured y');
        grid(ax,'on');
        xlabel(ax,'Time (s)');
        ylabel(ax,sprintf('h_%d (cm)',channel));
        title(ax,sprintf('%s MPC: h_%d',traceNames{col},channel));
        localStepLines(ax,config);
        legend(ax,'Location','best');
    end
end
for channel = 1:2
    limits = [inf -inf];
    for col = 1:2
        limits(1) = min(limits(1),axTrack(channel,col).YLim(1));
        limits(2) = max(limits(2),axTrack(channel,col).YLim(2));
    end
    for col = 1:2
        ylim(axTrack(channel,col),limits);
    end
end
figures = figTracking;

figDisturbance = experimentFigure('Name', ...
    'Problem 8: disturbance estimates','Position',[120 120 1000 650]);
sgtitle('Problem 8: F_3 and F_4 estimates with the biased MPC');
for channel = 1:2
    ax = subplot(2,1,channel);
    plot(ax,t,dTrue(channel,:),'k-','LineWidth',1.3, ...
        'DisplayName','true');
    hold(ax,'on');
    plot(ax,t,traceBiased.wHat(channel,:),'b--','LineWidth',1.3, ...
        'DisplayName','estimated');
    grid(ax,'on');
    xlabel(ax,'Time (s)');
    ylabel(ax,sprintf('\\DeltaF_%d (cm^3/s)',channel+2));
    title(ax,sprintf('Disturbance F_%d',channel+2));
    localStepLines(ax,config);
    legend(ax,'Location','best');
end
figures(end+1) = figDisturbance;

figPumps = experimentFigure('Name','Problem 8: physical pump flows', ...
    'Position',[140 140 1000 650]);
sgtitle('Problem 8: physical pump flows u_0 + u');
for channel = 1:2
    ax = subplot(2,1,channel);
    plot(ax,t,config.u0(channel) + traceUnbiased.u(channel,:),'b-', ...
        'LineWidth',1.2,'DisplayName','unbiased MPC');
    hold(ax,'on');
    plot(ax,t,config.u0(channel) + traceBiased.u(channel,:),'r-', ...
        'LineWidth',1.2,'DisplayName','biased MPC');
    grid(ax,'on');
    xlabel(ax,'Time (s)');
    ylabel(ax,sprintf('u_%d (cm^3/s)',channel));
    title(ax,sprintf('Pump %d',channel));
    localStepLines(ax,config);
    legend(ax,'Location','best');
end
figures(end+1) = figPumps;

figPrediction = experimentFigure('Name','Problem 8: horizon prediction', ...
    'Position',[160 160 1050 700]);
prediction = traceBiased.prediction;
predictionIndex = traceBiased.predictionIndex;
sgtitle(sprintf('Problem 8: biased MPC horizon prediction at t = %g s', ...
    t(predictionIndex)));
for channel = 1:2
    ax = subplot(2,1,channel);
    horizon = predictionIndex:min(predictionIndex + mpcBiased.N,Nsim);
    plot(ax,t(horizon),op.h(channel) + traceBiased.z(channel,horizon), ...
        'k-','LineWidth',1.4,'DisplayName','realized');
    hold(ax,'on');
    predictedSteps = 1:min(mpcBiased.N,Nsim - predictionIndex);
    plot(ax,t(predictionIndex + predictedSteps), ...
        op.h(channel) + prediction.Z(channel,predictedSteps), ...
        'b--o','LineWidth',1.3,'MarkerSize',4,'DisplayName','predicted');
    stairs(ax,t(horizon),op.h(channel) + r(channel,horizon), ...
        'm:','LineWidth',1.2,'DisplayName','reference');
    plot(ax,t(predictionIndex), ...
        op.h(channel) + traceBiased.z(channel,predictionIndex), ...
        'ko','MarkerFaceColor','k','HandleVisibility','off');
    xline(ax,t(predictionIndex),'k:','LineWidth',1.0, ...
        'HandleVisibility','off');
    grid(ax,'on');
    xlabel(ax,'Time (s)');
    ylabel(ax,sprintf('h_%d (cm)',channel));
    title(ax,sprintf('Tank %d level',channel));
    legend(ax,'Location','best');
end
figures(end+1) = figPrediction;

%% Summary
fprintf(['Problem 8 unbiased: RMSE [%.3f %.3f] cm, final offset ' ...
    '[%+.3f %+.3f] cm, max |e| %.3f cm\n'], ...
    metricsUnbiased.levelRMSE(1),metricsUnbiased.levelRMSE(2), ...
    metricsUnbiased.finalOffset(1),metricsUnbiased.finalOffset(2), ...
    max(metricsUnbiased.levelMaxAbs));
fprintf(['Problem 8 biased:   RMSE [%.3f %.3f] cm, final offset ' ...
    '[%+.3f %+.3f] cm, max |e| %.3f cm\n'], ...
    metricsBiased.levelRMSE(1),metricsBiased.levelRMSE(2), ...
    metricsBiased.finalOffset(1),metricsBiased.finalOffset(2), ...
    max(metricsBiased.levelMaxAbs));
fprintf(['Problem 8 effort:   unbiased |u|^2 = %.1f, |du|^2 = %.1f; ' ...
    'biased |u|^2 = %.1f, |du|^2 = %.1f\n'], ...
    metricsUnbiased.controlEffortTotal,metricsUnbiased.inputMoveEffortTotal, ...
    metricsBiased.controlEffortTotal,metricsBiased.inputMoveEffortTotal);
fprintf(['Problem 8 biased F estimates: final [%.2f %.2f], true ' ...
    '[%.2f %.2f] cm^3/s\n'], ...
    metricsBiased.disturbanceEstimateFinal(1), ...
    metricsBiased.disturbanceEstimateFinal(2), ...
    metricsBiased.disturbanceTrueFinal(1), ...
    metricsBiased.disturbanceTrueFinal(2));
if identifiedAvailable
    fprintf(['Problem 8 identified MPC: %d states, N = %d, ' ...
        'H positive definite = %d\n'], ...
        mpcIdentified.nx,mpcIdentified.N,mpcIdentified.Hpd);
else
    fprintf('Problem 8 identified MPC: model unavailable\n');
end
fprintf('Problem 8 prediction sample: t = %g s\n',t(predictionIndex));

%% Result
result = struct();
result.config = config;
result.model = models;
result.plant = struct('A',base.A,'B',base.B,'G',base.G,'C',base.C, ...
    'Ts',base.Ts,'op',op);
result.estimators = struct('unbiased',estUnbiased,'biased',estBiased);
result.mpc = struct('unbiased',mpcUnbiased,'biased',mpcBiased, ...
    'identified',mpcIdentified);
result.truth = struct('t',t,'d',dTrue,'w',w,'v',v,'ref',r);
result.controllers = struct('unbiased',traceUnbiased,'biased',traceBiased);
result.metrics = struct('unbiased',metricsUnbiased,'biased',metricsBiased);
result.prediction = struct('index',predictionIndex, ...
    'time',t(predictionIndex),'prediction',traceBiased.prediction);
end

function trace = localSimulate(mpc,base,dTrue,w,v,r,t,predictionIndex)
Nsim = numel(t);
trace = struct();
trace.t = t;
trace.x = zeros(mpc.nx,Nsim);
trace.z = zeros(mpc.nz,Nsim);
trace.y = zeros(mpc.nz,Nsim);
trace.u = zeros(mpc.nu,Nsim);
trace.xHat = zeros(mpc.nx,Nsim);
trace.wHat = zeros(mpc.nw,Nsim);
trace.innovation = zeros(mpc.nz,Nsim);
trace.prediction = [];
trace.predictionIndex = [];

x = base.x0;
xHat = mpc.estimator.x0;
if mpc.nw > 0
    wHat = mpc.estimator.w0;
else
    wHat = [];
end
uPrev = zeros(mpc.nu,1);

for k = 1:Nsim
    yk = base.C*x + v(:,k);
    zk = base.C*x;

    trace.x(:,k) = x;
    trace.z(:,k) = zk;
    trace.y(:,k) = yk;

    Rk = repmat(r(:,k),mpc.N,1);
    [uk,xHat,wHat,info] = mpc_step(mpc,xHat,wHat,uPrev,yk,Rk,[]);

    trace.u(:,k) = uk;
    trace.xHat(:,k) = xHat;
    trace.innovation(:,k) = info.innovation;
    if mpc.nw > 0
        trace.wHat(:,k) = wHat;
    end
    if k == predictionIndex
        trace.prediction = info.prediction;
        trace.predictionIndex = k;
    end

    if k < Nsim
        x = base.A*x + base.B*uk + base.G*(dTrue(:,k) + w(:,k));
    end
    uPrev = uk;
end
end

function m = localMetrics(trace,r,dTrue,config)
levelError = trace.z - r;
tracking = trace.t >= config.refTime;
if ~any(tracking)
    tracking = true(size(trace.t));
end
settle = trace.t >= trace.t(end) - config.settleTime;

m.levelError = levelError;
m.trackingIndex = tracking;
m.settleIndex = settle;
m.levelRMSE = sqrt(mean(levelError(:,tracking).^2,2));
m.levelMaxAbs = max(abs(levelError(:,tracking)),[],2);
m.finalOffset = mean(levelError(:,settle),2);
m.finalOffsetMaxAbs = max(abs(m.finalOffset));

inputMove = diff([zeros(size(trace.u,1),1), trace.u],1,2);
m.controlEffort = sum(trace.u.^2,2);
m.controlEffortTotal = sum(trace.u(:).^2);
m.inputMoveEffort = sum(inputMove.^2,2);
m.inputMoveEffortTotal = sum(inputMove(:).^2);
m.inputMin = min(trace.u,[],2);
m.inputMax = max(trace.u,[],2);
m.inputMoveMaxAbs = max(abs(inputMove),[],2);

if size(trace.wHat,1) > 0
    disturbanceError = trace.wHat - dTrue;
    m.disturbanceRMSE = sqrt(mean(disturbanceError.^2,2));
    m.disturbanceEstimateFinal = mean(trace.wHat(:,settle),2);
    m.disturbanceTrueFinal = dTrue(:,end);
    m.disturbanceErrorFinal = ...
        m.disturbanceEstimateFinal - m.disturbanceTrueFinal;
end
end

function localStepLines(ax,config)
yl = ylim(ax);
events = [config.refTime, config.stepTimes];
for k = 1:numel(events)
    plot(ax,[events(k) events(k)],yl,'k:','LineWidth',0.8, ...
        'HandleVisibility','off');
end
end

function s = localDefault(s,name,value)
if ~isfield(s,name) || isempty(s.(name))
    s.(name) = value;
end
end
