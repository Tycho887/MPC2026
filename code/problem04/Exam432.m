clear
clc
close all


dt = 1;            %tímaskilgreiningar  
Tend=20*60;
tgrid=0:dt:Tend;

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
F3 = zeros(1,N);
F4 = zeros(1,N);
F3(1)=20;    %arbitrary initial noise values 
F4(1)=20;

%% Disturbances for stochastic case
rng(67)


sigmaF3=4;
sigmaF4=4;

   
F3 = 20+sigmaF3*randn(1,N);
F4 = 20+sigmaF4*randn(1,N);

   

d=[F3;F4];


ds=[20;20];
us = [300;300];

xs = fsolve(@(x) ModifiedFourTankSystem(0,x,us,ds,p),xs_guess); %solving for steadystate 

x0=xs;

stepFactors = [1.10 1.25 1.50];


%storage
T_all = cell(2,3);
h1_norm_all = cell(2,3);
h2_norm_all = cell(2,3);
y1_norm_all = cell(2,3);
y2_norm_all = cell(2,3);


for inputCase = 1:2

for j = 1:3

F1 = 300*ones(1,N); %input fyrir F1 og F2 
F2 = 300*ones(1,N);

if inputCase == 1
    F1(tgrid(1:end-1)>=300)=300*stepFactors(j);  %step response at t=300
else
    F2(tgrid(1:end-1)>=300)=300*stepFactors(j);  %step response at t=300
end

U = [F1;F2];

x=x0;
T=tgrid(1);
X=x0.';

for k = 1:N

    dk = d(:,k);
    u = U(:,k); % all rows in column k

    [Tk,Xk] = ode45(@(t,x) ModifiedFourTankSystem(t,x,u,dk,p), ...
                    [tgrid(k) tgrid(k+1)], x);

    x = Xk(end,:)';
    T = [T;Tk(2:end);];
    X = [X;Xk(2:end,:);];

end


h1=X(:,1)/(rho*A1);
h2=X(:,2)/(rho*A2);
h3=X(:,3)/(rho*A3);
h4=X(:,4)/(rho*A4);

sigmaMeasurement = 0.10;

y1 = h1 + sigmaMeasurement*randn(size(h1));
y2 = h2 + sigmaMeasurement*randn(size(h2));

Z=[h1 h2];

%Normalizing
deltaF = 300*(stepFactors(j)-1);
h1_ss = xs(1)/(rho*A1);
h2_ss = xs(2)/(rho*A2);
y1_ss=h1_ss
y2_ss=h2_ss

h1_norm = (h1 - h1_ss)/deltaF;
h2_norm = (h2 - h2_ss)/deltaF;
y1_norm = (y1 - y1_ss)/deltaF;
y2_norm = (y2 - y2_ss)/deltaF;

%storing runs
T_all{inputCase,j} = T;
h1_norm_all{inputCase,j} = h1_norm;
h2_norm_all{inputCase,j} = h2_norm;
y1_norm_all{inputCase,j} = y1_norm;
y2_norm_all{inputCase,j} = y2_norm;



end

end


%%plot

figure

% F1 to y1
subplot(2,2,1)

plot(T_all{1,1},y1_norm_all{1,1},'LineWidth',1)
hold on
plot(T_all{1,2},y1_norm_all{1,2},'LineWidth',1)
plot(T_all{1,3},y1_norm_all{1,3},'LineWidth',1)

xlabel('Time [s]')
ylabel('\Delta y_1 / \Delta F_1')
title('F_1 to y_1')
legend('10%','25%','50%')
grid on


% F2 to y1
subplot(2,2,2)

plot(T_all{2,1},y1_norm_all{2,1},'LineWidth',1)
hold on
plot(T_all{2,2},y1_norm_all{2,2},'LineWidth',1)
plot(T_all{2,3},y1_norm_all{2,3},'LineWidth',1)

xlabel('Time [s]')
ylabel('\Delta y_1 / \Delta F_2')
title('F_2 to y_1')
legend('10%','25%','50%')
grid on


% F1 to y2
subplot(2,2,3)

plot(T_all{1,1},y2_norm_all{1,1},'LineWidth',1)
hold on
plot(T_all{1,2},y2_norm_all{1,2},'LineWidth',1)
plot(T_all{1,3},y2_norm_all{1,3},'LineWidth',1)

xlabel('Time [s]')
ylabel('\Delta y_2 / \Delta F_1')
title('F_1 to y_2')
legend('10%','25%','50%')
grid on


% F2 to y2
subplot(2,2,4)

plot(T_all{2,1},y2_norm_all{2,1},'LineWidth',1)
hold on
plot(T_all{2,2},y2_norm_all{2,2},'LineWidth',1)
plot(T_all{2,3},y2_norm_all{2,3},'LineWidth',1)

xlabel('Time [s]')
ylabel('\Delta y_2 / \Delta F_2')
title('F_2 to y_2')
legend('10%','25%','50%')
grid on


sgtitle('Medium noise')