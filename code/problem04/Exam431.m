clear
clc
close all

Ts = 60;            %tímaskilgreiningar  
Tend=20*60;
tgrid=0:Ts:Tend;

N = length(tgrid) -1;

F1 = 300*ones(1,N); %input fyrir F1 og F2 
F2 = 300*ones(1,N);

U = [F1;F2];

%% Parameters
a1 = 1.2272;
a2 = 1.2272;
a3 = 1.2272;
a4 = 1.2272;

A1 = 380.1327;
A2 = 380.1327;
A3 = 380.1327;
A4 = 380.1327;

gamma1 = 0.6;
gamma2 = 0.7;

g = 981;
rho = 1.00;

p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

%% Initial condition
xs_guess = [5000;
      5000;
      5000;
      5000];

%% Disturbances for deterministic case
F3 = 0;
F4 = 0;

d = [F3;
     F4];

us = [300;300];

xs = fsolve(@(x) ModifiedFourTankSystem(0,x,us,d,p),xs_guess); %solving for steadystate 

x0=xs;

stepFactors = [1.10 1.25 1.50];

T_all = cell(2,3);
h1_all = cell(2,3);
h2_all = cell(2,3);

for inputCase = 1:2

for j = 1:3

F1 = 300*ones(1,N); %input fyrir F1 og F2 
F2 = 300*ones(1,N);

if inputCase == 1
    F1(tgrid(1:end-1)>=300)=300*stepFactors(j);
else
    F2(tgrid(1:end-1)>=300)=300*stepFactors(j);
end;  
%step response at t=300

U = [F1;F2];

x=x0;
T=tgrid(1);
X=x0.';

for k = 1:N


    u = U(:,k); % all rows in column k

    [Tk,Xk] = ode45(@(t,x) ModifiedFourTankSystem(t,x,u,d,p),[tgrid(k) tgrid(k+1)], x);

    Tk(2:end)
    Xk(2:end,:)

    x = Xk(end,:)';
    T = [T;Tk(2:end)];
    X = [X; Xk(2:end,:)];

end

h1=X(:,1)/(rho*A1);
h2=X(:,2)/(rho*A2);
h3=X(:,3)/(rho*A3);
h4=X(:,4)/(rho*A4);

Z=[h1 h2];

%Normalizing
deltaF = 300*(stepFactors(j)-1);
h1_ss = xs(1)/(rho*A1);
h2_ss = xs(2)/(rho*A2);

h1_norm = (h1 - h1_ss)/deltaF;
h2_norm = (h2 - h2_ss)/deltaF;

%storing runs

T_all{inputCase,j} = T;
h1_all{inputCase,j} = h1_norm;
h2_all{inputCase,j} = h2_norm;

end
end



%%plot

figure

plot(T_all{1,1},h1_all{1,1},'LineWidth',1.5)
hold on
plot(T_all{1,2},h1_all{1,2},'LineWidth',1.5)
plot(T_all{1,3},h1_all{1,3},'LineWidth',1.5)

xlabel('Time [s]')
ylabel('\Delta h_1 / \Delta F_1')
legend('10%','25%','50%')
title('Normalized F_1 to h_1 step responses')
grid on



%% 4.4 Transfer function identification

% first index: input
% 1 = F1
% 2 = F2

% second index: step size
% 1 = 10%
% 2 = 25%
% 3 = 50%


stepTime = 300;
j = 1;                 % use 10% step

% F1 -> y1
t11 = T_all{1,j};
y11 = h1_all{1,j};

% F2 -> y1
t12 = T_all{2,j};
y12 = h1_all{2,j};

% F1 -> y2
t21 = T_all{1,j};
y21 = h2_all{1,j};

% F2 -> y2
t22 = T_all{2,j};
y22 = h2_all{2,j};


%using only parts after step
i11 = t11 >= stepTime;
tt11 = t11(i11) - stepTime;
yy11 = y11(i11);

i12 = t12 >= stepTime;
tt12 = t12(i12) - stepTime;
yy12 = y12(i12);

i21 = t21 >= stepTime;
tt21 = t21(i21) - stepTime;
yy21 = y21(i21);

i22 = t22 >= stepTime;
tt22 = t22(i22) - stepTime;
yy22 = y22(i22);

%theoretical
firstOrder = @(theta,t) theta(1)*(1-exp(-t/theta(2)));

%fitting f1 to y1
K11_guess = mean(yy11(end-20:end));

theta11 = lsqcurvefit(firstOrder, [K11_guess 70], tt11,yy11, [0 1],[Inf 1000]);

K11 = theta11(1);
T11 = theta11(2);


%fitting f2 to y2
K22_guess = mean(yy22(end-20:end));

theta22 = lsqcurvefit(firstOrder, [K22_guess 80], tt22,yy22, [0 1],[Inf 1000]);

K22 = theta22(1);
T22 = theta22(2);



secondOrder = @(theta,t) theta(1)*(1 -(theta(2)*exp(-t/theta(2)) - theta(3)*exp(-t/theta(3)))/(theta(2)-theta(3)));

%fitting F2  to y1
K12_guess = mean(yy12(end-20:end));

theta12 = lsqcurvefit(secondOrder, [K12_guess 70 25], tt12,yy12, [0 1 1],[Inf 1000 1000]);

K12 = theta12(1);
T12_1 = theta12(2);
T12_2 = theta12(3);

%fitting F1  to y2
K21_guess = mean(yy21(end-20:end));

theta21 = lsqcurvefit(secondOrder, [K21_guess 90 30], tt21,yy21, [0 1 1],[Inf 1000 1000]);

K21 = theta21(1);
T21_1 = theta21(2);
T21_2 = theta21(3);

fprintf('G11: K = %.5f, T = %.2f s\n',K11,T11)

fprintf('G12: K = %.5f, T1 = %.2f s, T2 = %.2f s\n', K12,T12_1,T12_2)

fprintf('G21: K = %.5f, T1 = %.2f s, T2 = %.2f s\n', K21,T21_1,T21_2)

fprintf('G22: K = %.5f, T = %.2f s\n',K22,T22)

%the transfer function

G11 = tf(K11,[T11 1]);

G12 = tf(K12,conv([T12_1 1],[T12_2 1]));

G21 = tf(K21,conv([T21_1 1],[T21_2 1]));

G22 = tf(K22,[T22 1]);

G = [G11 G12;
     G21 G22]

