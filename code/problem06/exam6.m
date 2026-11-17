
%% Problem 6 - Kalman filters
% Problem 5 ZOH and SDE models, and Problem 4 with both noise models
clear
clc
close all
rng(100)

%% Parameters and saved discrete time models
p = [1.2272*ones(4,1); 380.1327*ones(4,1); 0.6; 0.7; 981; 1];
gamma1 = p(9);
gamma2 = p(10);

load('problem5_statespace.mat', ...
    'Ad2','Bd2','Cd2','Ed2','Ad3','Bd3','Cd3','Qd3','Ts');
load('problem4_statespace.mat', ...
    'Ad_P4','Bd_P4','Cd_P4','Dd_P4','Ts_markov');

n4 = size(Ad_P4,1);
O4 = obsv(Ad_P4,Cd_P4);
H4 = zeros(2*n4,2);
Ak = eye(n4);

for j = 1:n4        %calculate markov parameters
    M = Cd_P4*Ak*Bd_P4;
    H4(2*j-1,1) = M(1,2)/(1-gamma2);
    H4(2*j,2)   = M(2,1)/(1-gamma1);
    Ak = Ak*Ad_P4;
end
E4 = O4\H4;

%% Noise covariances and uncertainties
% using same values as chap8
sigma_F = 5;                   % flow standard deviation
sigma_y = 0.5;                 % sensor standard deviation
Sigma_F = sigma_F^2*eye(2);    % noise covariance for F3 and F4
R = sigma_y^2*eye(2);          % sensor covariance


Q2 = Ed2*Sigma_F*Ed2';         % Problem 5 ZOH
Q3 = (Qd3+Qd3')/2;            % Problem 5 SDE
Q4 = E4*Sigma_F*E4';           % Problem 4 ZOH

% Problem 4, SDE
% sF as the SDE intensity
sF = 5;
sys4c = d2c(ss(Ad_P4,E4,Cd_P4,zeros(2),Ts_markov),'zoh');
Ac4 = sys4c.A;
Ec4 = sys4c.B;
W4 = Ec4*(sF^2*eye(2))*Ec4';
V4 = expm([-Ac4 W4; zeros(n4) Ac4']*Ts_markov);
Q4SDE = V4(n4+1:end,n4+1:end)'*V4(1:n4,n4+1:end);
Q4SDE = (Q4SDE+Q4SDE')/2;

P0 = 400^2*eye(4);             % Prior mass state uncertainty
P04 = eye(n4);                  % Prior realization state uncertainty
P0s = 10^2*eye(2);             % Prior unknown inflow uncertainty
Qs = 1^2*eye(2);               % Random walk covariance for inflows

%% Experiment definition
N = 500;
t5 = (0:N)*Ts;
t4 = (0:N)*Ts_markov;   %separate grid between problems

%pump inflow deviations and disturbances 
u5 = zeros(2,N);
u5(1,t5(1:N)>=150 & t5(1:N)<400) = 30;
u5(2,t5(1:N)>=450 & t5(1:N)<700) = 30;
u4 = zeros(2,N);
u4(1,t4(1:N)>=150 & t4(1:N)<400) = 30;
u4(2,t4(1:N)>=450 & t4(1:N)<700) = 30;

d5 = zeros(2,N);
d5(1,t5(1:N)>=300) = 15;
d5(2,t5(1:N)>=550) = -10;
d4 = zeros(2,N);
d4(1,t4(1:N)>=300) = 15;
d4(2,t4(1:N)>=550) = -10;


%random noises
x0 = [400;-350;200;-250];     
v5 = sigma_y*randn(2,N+1);
v4 = sigma_y*randn(2,N+1);
wF5 = sigma_F*randn(2,N);
wF4 = sigma_F*randn(2,N);
wSDE = chol(Q3,'lower')*randn(4,N);
%using eigenvalues, problem 4 SDE covariance can be singular
[V,D] = eig(Q4SDE);
w4SDE = V*diag(sqrt(max(diag(D),0)))*randn(n4,N);
D0 = zeros(2,2);

%% 6.2 Kalman filters without step changes
%separate measurements for each model
zero5 = zeros(2,N);
zero4 = zeros(2,N);

%Problem 5: ZOH 
[x2,y2] = simulateLinear(Ad2,Bd2,Cd2,D0,Ed2,u5,zero5,Ed2*wF5,v5,x0);
[xDyn2,xStat2] = kalmanPair(Ad2,Bd2,Cd2,D0,Q2,R,u5,y2,P0);

%Problem 5: SDE
[x3,y3] = simulateLinear(Ad3,Bd3,Cd3,D0,Ed2,u5,zero5,wSDE,v5,x0);
[xDyn3,xStat3] = kalmanPair(Ad3,Bd3,Cd3,D0,Q3,R,u5,y3,P0);

%Problem 4: ZOH 
[x4,y4] = simulateLinear(Ad_P4,Bd_P4,Cd_P4,Dd_P4,E4, u4,zero4,E4*wF4,v4,zeros(n4,1));
[xDyn4,xStat4] = kalmanPair(Ad_P4,Bd_P4,Cd_P4,Dd_P4, Q4,R,u4,y4,P04);

%Problem 4: SDE
[x4SDE,y4SDE] = simulateLinear(Ad_P4,Bd_P4,Cd_P4,Dd_P4,E4, u4,zero4,w4SDE,v4,zeros(n4,1));
[xDyn4SDE,xStat4SDE] = kalmanPair(Ad_P4,Bd_P4,Cd_P4,Dd_P4, Q4SDE,R,u4,y4SDE,P04);

fprintf('\n6.2: No step changes. Problem 5 level RMSE [cm]:\n')
disp(table((1:2)', sqrt(mean((Cd2*(xDyn2-x2)).^2,2)),sqrt(mean((Cd2*(xStat2-x2)).^2,2)), sqrt(mean((Cd3*(xDyn3-x3)).^2,2)),sqrt(mean((Cd3*(xStat3-x3)).^2,2)), 'VariableNames',{'Tank','P.5 (ZOH) Dynamic','P.5 (ZOH) Static','P.5 (SDE) Dynamic','P.5 (SDE) Static'}))

fprintf('\n6.2: No step changes. Problem 4 level RMSE [cm]:\n')
disp(table((1:2)', sqrt(mean((Cd_P4*(xDyn4-x4)).^2,2)), sqrt(mean((Cd_P4*(xStat4-x4)).^2,2)), sqrt(mean((Cd_P4*(xDyn4SDE-x4SDE)).^2,2)), sqrt(mean((Cd_P4*(xStat4SDE-x4SDE)).^2,2)), 'VariableNames',{'Tank','P.4 (ZOH) Dynamic','P.4 (ZOH) Static','P.4 (SDE) Dynamic','P.4 (SDE) Static'}))

plotLevelEstimates(t5,Cd2*x2,Cd2*xDyn2,Cd2*xStat2,'6.2 - Problem 5 (ZOH)');
plotLevelEstimates(t5,Cd3*x3,Cd3*xDyn3,Cd3*xStat3,'6.2 - Problem 5 (SDE)');
plotLevelEstimates(t4,Cd_P4*x4,Cd_P4*xDyn4,Cd_P4*xStat4, '6.2 - Problem 4 (ZOH)');
plotLevelEstimates(t4,Cd_P4*x4SDE,Cd_P4*xDyn4SDE,Cd_P4*xStat4SDE, '6.2 - Problem 4 (SDE)');

%% 6.3 Kalman filters with step changes
%Ordinary filters cannot explicitly estimate persistent unknown inflows
%In augmented models s(k+1)=s(k)+eta(k) with eta~N(0,Qs).

%Problem 5 ZOH: compare the usual filter with augmented one
[xStep2,yStep2] = simulateLinear(Ad2,Bd2,Cd2,D0,Ed2, u5,d5,Ed2*wF5,v5,x0);
[~,~,~,yBase2] = kalmanPair(Ad2,Bd2,Cd2,D0,Q2,R,u5,yStep2,P0);
[Aa2,Ba2,Ca2,Qa2,P0a2] = augmentInflows(Ad2,Bd2,Cd2,Ed2,Q2,Qs,P0,P0s);
[xAugDyn2,xAugStat2,yAugDyn2,yAugStat2] = kalmanPair(Aa2,Ba2,Ca2,D0,Qa2,R,u5,yStep2,P0a2);

%Problem 5 SDE: keep SDE covariance but estimate same physical inflow steps
[xStep3,yStep3] = simulateLinear(Ad3,Bd3,Cd3,D0,Ed2, u5,d5,wSDE,v5,x0);
[~,~,~,yBase3] = kalmanPair(Ad3,Bd3,Cd3,D0,Q3,R,u5,yStep3,P0);
[Aa3,Ba3,Ca3,Qa3,P0a3] = augmentInflows(Ad3,Bd3,Cd3,Ed2,Q3,Qs,P0,P0s);
[xAugDyn3,xAugStat3,yAugDyn3,yAugStat3] = kalmanPair(Aa3,Ba3,Ca3,D0,Qa3,R,u5,yStep3,P0a3);

% Problem 4 ZOH: augmented states estimate inflow steps
[xStep4,yStep4] = simulateLinear(Ad_P4,Bd_P4,Cd_P4,Dd_P4,E4, u4,d4,E4*wF4,v4,zeros(n4,1));
[~,~,~,yBase4] = kalmanPair(Ad_P4,Bd_P4,Cd_P4,Dd_P4, Q4,R,u4,yStep4,P04);
[Aa4,Ba4,Ca4,Qa4,P0a4] = augmentInflows(Ad_P4,Bd_P4,Cd_P4,E4,Q4,Qs,P04,P0s);
[xAugDyn4,xAugStat4,yAugDyn4,yAugStat4] = kalmanPair(Aa4,Ba4,Ca4,Dd_P4,Qa4,R,u4,yStep4,P0a4);

% Problem 4 SDE: same disturbance steps, but SDE process noise.
[xStep4SDE,yStep4SDE] = simulateLinear(Ad_P4,Bd_P4,Cd_P4,Dd_P4,E4, u4,d4,w4SDE,v4,zeros(n4,1));
[~,~,~,yBase4SDE] = kalmanPair(Ad_P4,Bd_P4,Cd_P4,Dd_P4, Q4SDE,R,u4,yStep4SDE,P04);
[Aa4SDE,Ba4SDE,Ca4SDE,Qa4SDE,P0a4SDE] = augmentInflows(Ad_P4,Bd_P4,Cd_P4,E4,Q4SDE,Qs,P04,P0s);
[xAugDyn4SDE,xAugStat4SDE,yAugDyn4SDE,yAugStat4SDE] = kalmanPair(Aa4SDE,Ba4SDE,Ca4SDE,Dd_P4,Qa4SDE,R,u4,yStep4SDE,P0a4SDE);

%prediction errors calculated BEFORE the current measurement update
%tail window starts long after both steps to assess remaining bias.
tail5 = find(t5 >= 0.8*t5(end));
tail4 = find(t4 >= 0.8*t4(end));
hStep2 = Cd2*xStep2;
hStep3 = Cd3*xStep3;
hStep4 = Cd_P4*xStep4 + Dd_P4*[u4 zeros(2,1)];
hStep4SDE = Cd_P4*xStep4SDE + Dd_P4*[u4 zeros(2,1)];

fprintf('\n6.3: Average one step level prediction error [cm], Problem 5 ZOH:\n')
showBias(hStep2,yBase2,yAugDyn2,yAugStat2,tail5)
fprintf('\n6.3: Average one step level prediction error [cm], Problem 5 SDE:\n')
showBias(hStep3,yBase3,yAugDyn3,yAugStat3,tail5)
fprintf('\n6.3: Average one step level prediction error [cm], Problem 4 ZOH:\n')
showBias(hStep4,yBase4,yAugDyn4,yAugStat4,tail4)
fprintf('\n6.3: Average one step level prediction error [cm], Problem 4 SDE:\n')
showBias(hStep4SDE,yBase4SDE,yAugDyn4SDE,yAugStat4SDE,tail4)

plotPredictions(t5,hStep2,yBase2,yAugDyn2,yAugStat2, '6.3 - Problem 5 ZOH : one step output predictions');
plotPredictions(t5,hStep3,yBase3,yAugDyn3,yAugStat3, '6.3 - Problem 5 SDE: one step output predictions');
plotPredictions(t4,hStep4,yBase4,yAugDyn4,yAugStat4, '6.3 - Problem 4 ZOH: one step output predictions');
plotPredictions(t4,hStep4SDE,yBase4SDE,yAugDyn4SDE,yAugStat4SDE, '6.3 - Problem 4 SDE: one step output predictions');

figure
for i = 1:2
    subplot(2,1,i)
    stairs(t5(1:N),d5(i,:),'k','LineWidth',1.4)
    hold on
    plot(t5,xAugDyn2(4+i,:),'b','LineWidth',1.1)
    plot(t5,xAugStat2(4+i,:),'r--','LineWidth',1.1)
    xlabel('Time [s]')
    ylabel(['\Delta F_' num2str(i+2) ' [cm^3/s]'])
    legend('True step','Augmented dynamic','Augmented static')
    grid on
end
sgtitle('6.3 - Problem 5 ZOH: estimated unmeasured inflows')

%% 6.4 Model 2 filters on the nonlinear plant
% Use physical nonlinear plant, but keep same linear filter
us = [300;300];
ds = [20;20];
xs = fsolve(@(x) ModifiedFourTankSystem(0,x,us,ds,p),5000*ones(4,1));

xNL = zeros(4,N+1);
xNL(:,1) = xs+x0;
for k = 1:N
    uk = us+u5(:,k);
    dk = ds+wF5(:,k);
    [~,X] = ode45(@(t,x) ModifiedFourTankSystem(t,x,uk,dk,p), [t5(k) t5(k+1)],xNL(:,k));
    xNL(:,k+1) = X(end,:)';
end
xDevNL = xNL-xs;
yNL = Cd2*xDevNL+v5;
[xDynNL,xStatNL] = kalmanPair(Ad2,Bd2,Cd2,D0,Q2,R,u5,yNL,P0);

fprintf('\n6.4: Nonlinear plant, Model 2 level RMSE [cm]:\n')
disp(table((1:2)',sqrt(mean((Cd2*(xDynNL-xDevNL)).^2,2)), sqrt(mean((Cd2*(xStatNL-xDevNL)).^2,2)), 'VariableNames',{'Tank','Dynamic','Static'}))
plotLevelEstimates(t5,Cd2*xDevNL,Cd2*xDynNL,Cd2*xStatNL, '6.4 - Problem 5 ZOH Kalman filters on the nonlinear plant');

%% Functions
function [x,y] = simulateLinear(A,B,C,D,E,u,d,w,v,x0)
    N = size(u,2);
    x = zeros(size(A,1),N+1);
    y = zeros(size(C,1),N+1);
    x(:,1) = x0;
    for k = 1:N
        y(:,k) = C*x(:,k)+D*u(:,k)+v(:,k);
        x(:,k+1) = A*x(:,k)+B*u(:,k)+E*d(:,k)+w(:,k);
    end
    y(:,N+1) = C*x(:,N+1)+v(:,N+1);
end

function [xDyn,xStat,yPredDyn,yPredStat] = kalmanPair(A,B,C,D,Q,R,u,y,P0)
    N = size(u,2);
    nx = size(A,1);
    ny = size(C,1);
    xDyn = zeros(nx,N+1);
    xStat = zeros(nx,N+1);
    yPredDyn = zeros(ny,N+1);
    yPredStat = zeros(ny,N+1);

    Pss = dare(A',C',Q,R);                %steady state prior covariance
    Kstat = (Pss*C')/(C*Pss*C'+R);       %constant measurement update gain
    Pprior = P0;
    xPriorDyn = zeros(nx,1);
    xPriorStat = zeros(nx,1);

    for k = 1:N+1
        if k <= N
            uk = u(:,k);
        else
            uk = zeros(size(B,2),1);
        end
        yPredDyn(:,k) = C*xPriorDyn+D*uk;
        yPredStat(:,k) = C*xPriorStat+D*uk;

        Re = C*Pprior*C'+R;
        Kdyn = (Pprior*C')/Re;
        xDyn(:,k) = xPriorDyn+Kdyn*(y(:,k)-yPredDyn(:,k));
        xStat(:,k) = xPriorStat+Kstat*(y(:,k)-yPredStat(:,k));
        Ppost = Pprior-Kdyn*Re*Kdyn';

        if k <= N
            xPriorDyn = A*xDyn(:,k)+B*uk;
            xPriorStat = A*xStat(:,k)+B*uk;
            Pprior = A*Ppost*A'+Q;  
        end
    end
end

function [Aa,Ba,Ca,Qa,P0a] = augmentInflows(A,B,C,E,Q,Qs,P0,P0s)
    nx = size(A,1);
    nd = size(E,2);
    Aa = [A E; zeros(nd,nx) eye(nd)];
    Ba = [B; zeros(nd,size(B,2))];
    Ca = [C zeros(size(C,1),nd)];
    Qa = blkdiag(Q,Qs);
    P0a = blkdiag(P0,P0s);
end

function showBias(h,yBase,yDyn,yStat,tail)
    disp(table((1:2)',mean(h(:,tail)-yBase(:,tail),2), ...
        mean(h(:,tail)-yDyn(:,tail),2), ...
        mean(h(:,tail)-yStat(:,tail),2), ...
        'VariableNames',{'Tank','OrdinaryStatic','AugDynamic','AugStatic'}))
end

function plotStateEstimates(t,x,xDyn,xStat,heading)
    figure
    for i = 1:size(x,1)
        subplot(2,2,i)
        plot(t,x(i,:),'k','LineWidth',1.3)
        hold on
        plot(t,xDyn(i,:),'b','LineWidth',1.1)
        plot(t,xStat(i,:),'r--','LineWidth',1.1)
        xlabel('Time [s]')
        ylabel(['\Delta m_' num2str(i) ' [g]'])
        legend('True','Dynamic','Static')
        grid on
    end
    sgtitle(heading)
end

function plotLevelEstimates(t,h,hDyn,hStat,heading)
    figure
    for i = 1:2
        subplot(2,1,i)
        plot(t,h(i,:),'k','LineWidth',1.3)
        hold on
        plot(t,hDyn(i,:),'b','LineWidth',1.1)
        plot(t,hStat(i,:),'r--','LineWidth',1.1)
        xlabel('Time [s]')
        ylabel(['\Delta h_' num2str(i) ' [cm]'])
        legend('True','Dynamic','Static')
        grid on
    end
    sgtitle(heading)
end

function plotPredictions(t,h,yBase,yDyn,yStat,heading)
    figure
    for i = 1:2
        subplot(2,1,i)
        plot(t,h(i,:),'k','LineWidth',1.3)
        hold on
        plot(t,yBase(i,:),'r--','LineWidth',1.1)
        plot(t,yDyn(i,:),'b','LineWidth',1.1)
        plot(t,yStat(i,:),'g--','LineWidth',1.1)
        xlabel('Time [s]')
        ylabel(['\Delta h_' num2str(i) ' [cm]'])
        legend('True','Ordinary static','Augmented dynamic','Augmented static')
        grid on
    end
    sgtitle(heading)
end
