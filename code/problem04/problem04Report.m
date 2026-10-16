function [result, figures] = problem04Report(deterministic, stochastic, identified, markov)
%PROBLEM04REPORT Normalize all step tests and quantify report model accuracy.
% Uses existing trajectories; does not draw new noise or refit the models.
figures = gobjects(0); colors = lines(3);
result = struct('deterministicGains',zeros(2,2,3), ...
    'validationRMSE',zeros(2,2,3),'fitRMSE',zeros(2), ...
    'identifiedGains',dcgain(identified.G),'backend',identified.backend);
for j = 1:2
    figures(end+1) = experimentFigure('Name',sprintf('Normalized deterministic F%d',j),'Color','w');
    for a = 1:3
        r = deterministic.experiments{j,a};
        amplitude = deterministic.op.u(j)*deterministic.fractions(a);
        response = (r.z-deterministic.op.z)/amplitude;
        for i = 1:2
            subplot(2,1,i); hold on; grid on;
            plot(r.t-deterministic.stepTime,response(i,:), ...
                'Color',colors(a,:),'LineWidth',1.4,'DisplayName',sprintf('+%g%%',100*deterministic.fractions(a)));
            ylabel(sprintf('S_{%d%d} (cm/(cm^3/s))',i,j));
            xlabel('Time since step (s)'); xlim([-200 1000]); legend('Location','best');
            result.deterministicGains(i,j,a) = response(i,end);
            % Independent deterministic validation on a uniform 1 s grid.
            tv = (0:1000)';
            truth = interp1(r.t-deterministic.stepTime,response(i,:),tv,'pchip');
            prediction = step(identified.channels{i,j},tv);
            result.validationRMSE(i,j,a) = sqrt(mean((truth-prediction).^2));
        end
    end
end
for j = 1:2
    figures(end+1) = experimentFigure('Name',sprintf('Normalized stochastic F%d',j), ...
        'Color','w','Position',[100 100 1100 850]);
    for b = 1:3
        for a = 1:3
            r = stochastic.experiments{j,b,a};
            amplitude = stochastic.op.u(j)*stochastic.fractions(a);
            pre = r.t>=150 & r.t<200;
            baseline = mean(r.y(:,pre),2);
            for i = 1:2
                subplot(3,2,2*(b-1)+i); hold on; grid on;
                plot(r.t-200,(r.z(i,:)-stochastic.op.z(i))/amplitude, ...
                    'Color',colors(a,:),'LineWidth',1.2,'DisplayName',sprintf('+%g%%',100*stochastic.fractions(a)));
                idx = 1:10:numel(r.t);
                plot(r.t(idx)-200,(r.y(i,idx)-baseline(i))/amplitude,'.', ...
                    'Color',colors(a,:),'MarkerSize',3,'HandleVisibility','off');
                title(sprintf('h%d, noise scale %g',i,stochastic.noiseScales(b)));
                ylabel(sprintf('S_{%d%d} (cm/(cm^3/s))',i,j));
                xlabel('Time since step (s)'); xlim([-200 1000]);
                if b==1, legend('Location','best'); end
            end
        end
    end
end
for i = 1:2
    for j = 1:2
        n = identified.normalized{i,j}; mask = n.t>=identified.stepTime;
        result.fitRMSE(i,j) = sqrt(mean((n.true(mask)-n.model(mask)).^2));
    end
end
result.markovFirst = markov.H(:,:,2);
result.markovPartialSum = sum(markov.H,3);
result.markovTailFraction = 1-result.markovPartialSum./result.identifiedGains;
end
