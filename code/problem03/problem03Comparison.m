function [result, figures] = problem03Comparison()
% Independent mass-feedback closed loops with identical measurement noise.
p = fourTankParameters(); op = fourTankSteadyState([300;300],[0;0],p);
model = linearizeFourTank(op,p); [Gm,~] = fourTankTransferFunctions(model);
controller = tuneFourTankPI(Gm);
dt = 0.5; t = (0:dt:2000)'; N = numel(t);
reference = repmat(op.x(1:2),1,N);
reference(1,t>=200) = reference(1,t>=200)+2000;
reference(2,t>=1000) = reference(2,t>=1000)+2000;
noiseLevels = [0 10 50]; seed = 100;
rng(seed,'twister'); standardNoise = randn(2,N);
kp = [controller.C1.Kp;controller.C2.Kp];
ki = [controller.C1.Ki;controller.C2.Ki];
experiments = cell(1,numel(noiseLevels));
figures = experimentFigure('Name','Independent linear and nonlinear closed loops','Color','w');
colors = lines(numel(noiseLevels));
for j = 1:numel(noiseLevels)
    noise = noiseLevels(j)*standardNoise;
    xNL = repmat(op.x,1,N); xLinear = zeros(4,N);
    integralNL = zeros(2,1); integralLinear = zeros(2,1);
    uNL = zeros(2,N); uLinear = zeros(2,N);
    for k = 1:N
        errorNL = reference(:,k)-xNL(1:2,k)-noise(:,k);
        errorLinear = reference(:,k)-op.x(1:2)-xLinear(1:2,k)-noise(:,k);
        integralNL = integralNL + dt*errorNL;
        integralLinear = integralLinear + dt*errorLinear;
        uNL(:,k) = min(max(op.u+controller.W*(kp.*errorNL+ki.*integralNL),0),1000);
        uLinear(:,k) = min(max(op.u+controller.W*(kp.*errorLinear+ki.*integralLinear),0),1000);
        if k < N
            xNL(:,k+1) = xNL(:,k)+dt*fourTankDynamics(t(k),xNL(:,k),uNL(:,k),op.d,p);
            xLinear(:,k+1) = xLinear(:,k)+dt*(model.A*xLinear(:,k)+model.B*(uLinear(:,k)-op.u));
        end
    end
    xLinearAbsolute = xLinear+op.x;
    levelNoise = noise ./ (p.rho*p.A(1:2));
    experiments{j} = struct('noiseStdGrams',noiseLevels(j),'noiseGrams',noise, ...
        'xNonlinear',xNL,'xLinearDeviation',xLinear,'xLinear',xLinearAbsolute, ...
        'uNonlinear',uNL,'uLinear',uLinear, ...
        'zNonlinear',fourTankOutput(xNL,p),'zLinear',fourTankOutput(xLinearAbsolute,p), ...
        'yNonlinear',fourTankSensor(xNL,p,levelNoise), ...
        'yLinear',fourTankSensor(xLinearAbsolute,p,levelNoise));
    for tank = 1:2
        subplot(2,1,tank); hold on;
        plot(t,xNL(tank,:),'Color',colors(j,:),'DisplayName',sprintf('Nonlinear, sigma=%g g',noiseLevels(j)));
        plot(t,xLinearAbsolute(tank,:),'--','Color',colors(j,:),'DisplayName',sprintf('Linear, sigma=%g g',noiseLevels(j)));
    end
end
for tank = 1:2
    subplot(2,1,tank); plot(t,reference(tank,:),'k:','DisplayName','Reference');
    grid on; xlabel('Time (s)'); ylabel('Mass (g)'); title(sprintf('Tank %d mass tracking',tank)); legend('Location','best');
end
result = struct('p',p,'op',op,'model',model,'controller',controller,'t',t, ...
    'referenceGrams',reference,'seed',seed,'sampleTime',dt,'solver','euler', ...
    'pumpBounds',[0;1000],'experiments',{experiments});
end
