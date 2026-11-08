function [result, figures] = problem06Simulate(caseName, config)
%PROBLEM06SIMULATE Simulate and evaluate the Problem 6.2/6.3 filters.
%
%   [result, figures] = problem06Simulate('6.2')
%   [result, figures] = problem06Simulate('6.3')
%
% Case 6.2 uses the 4-state mass model and zero-mean stochastic disturbances.
% Case 6.3 uses the 6-state disturbance-augmented model and disturbances with
% step changes. Both cases run a static and a dynamic Kalman filter on the
% same truth trajectory and measurements.

if nargin < 1 || isempty(caseName)
    caseName = '6.2';
end
if nargin < 2 || isempty(config)
    config = struct();
end
caseName = validatestring(caseName,{'6.2','6.3'},mfilename,'caseName');

config = localDefault(config,'seed',6002 + strcmp(caseName,'6.3'));
config = localDefault(config,'tFinal',1200);
config = localDefault(config,'stepTimes',[300 600]);
config = localDefault(config,'stepSizes',[20 20]);
config = localDefault(config,'settleTime',200);

model = problem06Model(config);
base = model.base;
aug = model.augmented;
isStepCase = strcmp(caseName,'6.3');

%% Truth trajectory and measurements
Ts = base.Ts;
t = (0:Ts:config.tFinal)';
N = numel(t);
u = zeros(2,N);                         % pump-flow deviations
dTrue = zeros(2,N);                     % disturbance bias (deviation)
if isStepCase
    dTrue(1,t >= config.stepTimes(1)) = config.stepSizes(1);
    dTrue(2,t >= config.stepTimes(2)) = config.stepSizes(2);
end

rng(config.seed,'twister');
w = model.config.sigmaF * randn(2,N);
v = model.config.sigmaV * randn(2,N);

x = zeros(4,N);
x(:,1) = base.x0;
for k = 1:N-1
    x(:,k+1) = base.A*x(:,k) + base.B*u(:,k) ...
        + base.G*(dTrue(:,k) + w(:,k));
end
y = base.C*x + v;
z = base.C*x;

truth = struct('t',t,'u',u,'d',dTrue,'w',w,'v',v,'x',x,'y',y,'z',z);

%% Filters
[estStaticBase,filterStaticBase] = problem06Filter(base,'static',y,u);
[estDynamicBase,filterDynamicBase] = problem06Filter(base,'dynamic',y,u);

if isStepCase
    [estStaticAug,filterStaticAug] = problem06Filter(aug,'static',y,u);
    [estDynamicAug,filterDynamicAug] = problem06Filter(aug,'dynamic',y,u);
    estimates = struct('static',estStaticAug,'dynamic',estDynamicAug, ...
        'baselineStatic',estStaticBase);
    filters = struct('static',filterStaticAug,'dynamic',filterDynamicAug, ...
        'baselineStatic',filterStaticBase);
else
    estimates = struct('static',estStaticBase,'dynamic',estDynamicBase);
    filters = struct('static',filterStaticBase,'dynamic',filterDynamicBase);
end

%% Metrics
metrics = struct();
if isStepCase
    metrics.static = localMetrics(estStaticAug,z,x,dTrue,t,config,true,true);
    metrics.dynamic = localMetrics(estDynamicAug,z,x,dTrue,t,config,true,true);
    metrics.baselineStatic = localMetrics(estStaticBase,z,x,dTrue,t,config,true,false);
else
    metrics.static = localMetrics(estStaticBase,z,x,[],t,config,false,false);
    metrics.dynamic = localMetrics(estDynamicBase,z,x,[],t,config,false,false);
end

%% Figures
figures = gobjects(0);
hTrue = model.op.h(1:2) + z;
hMeas = model.op.h(1:2) + y;
hStaticBase = model.op.h(1:2) + estStaticBase.zFilt;
hDynamicBase = model.op.h(1:2) + estDynamicBase.zFilt;

if isStepCase
    hStaticAug = model.op.h(1:2) + estStaticAug.zFilt;
    hDynamicAug = model.op.h(1:2) + estDynamicAug.zFilt;

    fig = experimentFigure('Name','Problem 6.3: Kalman filter estimates', ...
        'Position',[100 100 1100 1000]);
    sgtitle('Problem 6.3: stochastic disturbances with step changes');

    localPlotPair(subplot(4,2,1),t,hTrue(1,:),hStaticAug(1,:),hDynamicAug(1,:), ...
        hMeas(1,:),'h_1 (cm)','Tank 1 level',config.stepTimes);
    localPlotPair(subplot(4,2,2),t,hTrue(2,:),hStaticAug(2,:),hDynamicAug(2,:), ...
        hMeas(2,:),'h_2 (cm)','Tank 2 level',config.stepTimes);

    localPlotPair(subplot(4,2,3),t,x(1,:),estStaticAug.xFilt(1,:), ...
        estDynamicAug.xFilt(1,:),[],'\Delta m_1 (g)','Tank 1 mass',config.stepTimes);
    localPlotPair(subplot(4,2,4),t,x(2,:),estStaticAug.xFilt(2,:), ...
        estDynamicAug.xFilt(2,:),[],'\Delta m_2 (g)','Tank 2 mass',config.stepTimes);
    localPlotPair(subplot(4,2,5),t,x(3,:),estStaticAug.xFilt(3,:), ...
        estDynamicAug.xFilt(3,:),[],'\Delta m_3 (g)','Tank 3 mass',config.stepTimes);
    localPlotPair(subplot(4,2,6),t,x(4,:),estStaticAug.xFilt(4,:), ...
        estDynamicAug.xFilt(4,:),[],'\Delta m_4 (g)','Tank 4 mass',config.stepTimes);

    localPlotPair(subplot(4,2,7),t,dTrue(1,:),estStaticAug.xFilt(5,:), ...
        estDynamicAug.xFilt(5,:),[],'\Delta F_3 (cm^3/s)', ...
        'Disturbance F_3',config.stepTimes);
    localPlotPair(subplot(4,2,8),t,dTrue(2,:),estStaticAug.xFilt(6,:), ...
        estDynamicAug.xFilt(6,:),[],'\Delta F_4 (cm^3/s)', ...
        'Disturbance F_4',config.stepTimes);
    figures(end+1) = fig;

    % The 4-state filter has no disturbance state, so it cannot remove the
    % offset caused by the step. This figure makes that comparison explicit.
    figBaseline = experimentFigure('Name','Problem 6.3: 4-state baseline offset', ...
        'Position',[150 150 1000 520]);
    sgtitle('Problem 6.3: 4-state filter (no disturbance state) vs 6-state filter');
    localPlotError(subplot(2,1,1),t,estStaticBase.zFilt(1,:)-z(1,:), ...
        estStaticAug.zFilt(1,:)-z(1,:),'h_1 error (cm)','Tank 1 level error',config.stepTimes);
    localPlotError(subplot(2,1,2),t,estStaticBase.zFilt(2,:)-z(2,:), ...
        estStaticAug.zFilt(2,:)-z(2,:),'h_2 error (cm)','Tank 2 level error',config.stepTimes);
    figures(end+1) = figBaseline;
else
    fig = experimentFigure('Name','Problem 6.2: Kalman filter estimates', ...
        'Position',[100 100 1100 850]);
    sgtitle('Problem 6.2: stochastic disturbances without step changes');

    localPlotPair(subplot(3,2,1),t,hTrue(1,:),hStaticBase(1,:),hDynamicBase(1,:), ...
        hMeas(1,:),'h_1 (cm)','Tank 1 level',[]);
    localPlotPair(subplot(3,2,2),t,hTrue(2,:),hStaticBase(2,:),hDynamicBase(2,:), ...
        hMeas(2,:),'h_2 (cm)','Tank 2 level',[]);
    localPlotPair(subplot(3,2,3),t,x(1,:),estStaticBase.xFilt(1,:), ...
        estDynamicBase.xFilt(1,:),[],'\Delta m_1 (g)','Tank 1 mass',[]);
    localPlotPair(subplot(3,2,4),t,x(2,:),estStaticBase.xFilt(2,:), ...
        estDynamicBase.xFilt(2,:),[],'\Delta m_2 (g)','Tank 2 mass',[]);
    localPlotPair(subplot(3,2,5),t,x(3,:),estStaticBase.xFilt(3,:), ...
        estDynamicBase.xFilt(3,:),[],'\Delta m_3 (g)','Tank 3 mass',[]);
    localPlotPair(subplot(3,2,6),t,x(4,:),estStaticBase.xFilt(4,:), ...
        estDynamicBase.xFilt(4,:),[],'\Delta m_4 (g)','Tank 4 mass',[]);
    figures(end+1) = fig;
end

%% Summary
if isStepCase
    fprintf(['Problem 6.3 static:  level RMSE [%.3f %.3f] cm, ' ...
        'offset max %.3f cm\n'], ...
        metrics.static.levelRMSE(1),metrics.static.levelRMSE(2), ...
        metrics.static.levelOffsetMaxAbs);
    fprintf(['Problem 6.3 dynamic: level RMSE [%.3f %.3f] cm, ' ...
        'offset max %.3f cm\n'], ...
        metrics.dynamic.levelRMSE(1),metrics.dynamic.levelRMSE(2), ...
        metrics.dynamic.levelOffsetMaxAbs);
    fprintf(['Problem 6.3 baseline 4-state: level RMSE [%.3f %.3f] cm, ' ...
        'offset max %.3f cm\n'], ...
        metrics.baselineStatic.levelRMSE(1), ...
        metrics.baselineStatic.levelRMSE(2), ...
        metrics.baselineStatic.levelOffsetMaxAbs);
else
    fprintf(['Problem 6.2 static:  level RMSE [%.3f %.3f] cm, ' ...
        'offset max %.3f cm\n'], ...
        metrics.static.levelRMSE(1),metrics.static.levelRMSE(2), ...
        metrics.static.levelOffsetMaxAbs);
    fprintf(['Problem 6.2 dynamic: level RMSE [%.3f %.3f] cm, ' ...
        'offset max %.3f cm\n'], ...
        metrics.dynamic.levelRMSE(1),metrics.dynamic.levelRMSE(2), ...
        metrics.dynamic.levelOffsetMaxAbs);
end

result = struct();
result.caseName = caseName;
result.config = config;
result.model = model;
result.truth = truth;
result.estimates = estimates;
result.filters = filters;
result.metrics = metrics;
end

function m = localMetrics(est,z,x,dTrue,t,config,hasSteps,hasDisturbanceState)
levelError = est.zFilt - z;
massError = est.xFilt(1:4,:) - x;

m.levelRMSE = sqrt(mean(levelError.^2,2));
m.massRMSE = sqrt(mean(massError.^2,2));
m.levelMaxAbs = max(abs(levelError),[],2);
m.massMaxAbs = max(abs(massError),[],2);

if hasSteps
    idx1 = t >= config.stepTimes(1) + config.settleTime & ...
        t < config.stepTimes(2);
    idx2 = t >= config.stepTimes(2) + config.settleTime;
    m.levelOffset = [mean(levelError(:,idx1),2), mean(levelError(:,idx2),2)];
    m.massOffset = [mean(massError(:,idx1),2), mean(massError(:,idx2),2)];
else
    idx = t >= t(end) - config.settleTime;
    m.levelOffset = mean(levelError(:,idx),2);
    m.massOffset = mean(massError(:,idx),2);
end
m.levelOffsetMaxAbs = max(abs(m.levelOffset),[],'all');
m.massOffsetMaxAbs = max(abs(m.massOffset),[],'all');

if hasDisturbanceState
    dError = est.xFilt(5:6,:) - dTrue;
    m.disturbanceRMSE = sqrt(mean(dError.^2,2));
    m.disturbanceMaxAbs = max(abs(dError),[],2);
    if hasSteps
        m.disturbanceOffset = [mean(dError(:,idx1),2), mean(dError(:,idx2),2)];
    end
end
end

function localPlotPair(ax,t,truth,staticEst,dynamicEst,meas,ylabelText,titleText,stepTimes)
plot(ax,t,truth,'k-','LineWidth',1.3,'DisplayName','true');
hold(ax,'on');
plot(ax,t,staticEst,'b--','LineWidth',1.2,'DisplayName','static KF');
plot(ax,t,dynamicEst,'r-.','LineWidth',1.2,'DisplayName','dynamic KF');
if ~isempty(meas)
    stride = max(1,round(numel(t)/50));
    idx = 1:stride:numel(t);
    plot(ax,t(idx),meas(idx),'.','Color',[0.45 0.45 0.45], ...
        'MarkerSize',7,'HandleVisibility','off');
end
grid(ax,'on');
xlabel(ax,'Time (s)');
ylabel(ax,ylabelText);
title(ax,titleText);
localStepLines(ax,t,stepTimes);
legend(ax,'Location','best');
end

function localPlotError(ax,t,baselineError,augmentedError,ylabelText,titleText,stepTimes)
plot(ax,t,baselineError,'m-','LineWidth',1.2,'DisplayName','4-state static KF');
hold(ax,'on');
plot(ax,t,augmentedError,'b--','LineWidth',1.2,'DisplayName','6-state static KF');
grid(ax,'on');
xlabel(ax,'Time (s)');
ylabel(ax,ylabelText);
title(ax,titleText);
localStepLines(ax,t,stepTimes);
legend(ax,'Location','best');
end

function localStepLines(ax,~,stepTimes)
if isempty(stepTimes)
    return;
end
yl = ylim(ax);
for k = 1:numel(stepTimes)
    plot(ax,[stepTimes(k) stepTimes(k)],yl,'k:','LineWidth',0.8, ...
        'HandleVisibility','off');
end
end

function s = localDefault(s,name,value)
if ~isfield(s,name) || isempty(s.(name))
    s.(name) = value;
end
end
