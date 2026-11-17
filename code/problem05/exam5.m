clear
clc
close all

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

%% 5.1 Linearization - Model 1

us = [300;300];
ds = [0;0];

%% Initial condition
xs_guess = [5000;
      5000;
      5000;
      5000];


xs = fsolve(@(x) ModifiedFourTankSystem(0,x,us,ds,p),xs_guess);

%steady-state heights
At = [A1;A2;A3;A4];
ap = [a1;a2;a3;a4];

hs = xs./(rho*At);

%tank time constants
Ttank = (At./ap).*sqrt(2*hs/g);

A_lin1 = [-1/Ttank(1)  0             1/Ttank(3)   0;
           0            -1/Ttank(2)  0             1/Ttank(4);
           0             0           -1/Ttank(3)   0;
           0             0            0            -1/Ttank(4)];

B_lin1 = [rho*gamma1      0;
          0               rho*gamma2;
          0               rho*(1-gamma2);
          rho*(1-gamma1)  0];

E_lin1 = [0    0;
          0    0;
          rho  0;
          0    rho];

C_lin1 = [1/(rho*A1)  0            0  0;
          0           1/(rho*A2)   0  0];

D_lin1 = zeros(2,2);

A_lin1
B_lin1
E_lin1
C_lin1
D_lin1

%% 5.1 Model 2 - piecewise constant stochastic disturbances

us2 = [300;300];
ds2 = [20;20];

xs2 = fsolve(@(x) ModifiedFourTankSystem(0,x,us2,ds2,p),xs_guess);

%steady-state heights
hs2 = xs2./(rho*At);

%tank time constants
Ttank2 = (At./ap).*sqrt(2*hs2/g);

A_lin2 = [-1/Ttank2(1)  0               1/Ttank2(3)   0;
           0             -1/Ttank2(2)   0              1/Ttank2(4);
           0              0             -1/Ttank2(3)   0;
           0              0              0             -1/Ttank2(4)];

B_lin2 = [rho*gamma1      0;
          0               rho*gamma2;
          0               rho*(1-gamma2);
          rho*(1-gamma1)  0];

E_lin2 = [0    0;
          0    0;
          rho  0;
          0    rho];

C_lin2 = [1/(rho*A1)  0            0  0;
          0           1/(rho*A2)   0  0];

D_lin2 = zeros(2,2);

xs2
hs2
Ttank2
A_lin2
B_lin2
E_lin2
C_lin2
D_lin2

%% 5.2 Continuous-time transfer functions

%Model 1
sysc1 = ss(A_lin1,B_lin1,C_lin1,D_lin1);
Glin1 = minreal(tf(sysc1));

disp('Model 1 transfer functions:')
Glin1

%Model 2
sysc2 = ss(A_lin2,B_lin2,C_lin2,D_lin2);
Glin2 = minreal(tf(sysc2));

disp('Model 2 transfer functions:')
Glin2

%% 5.1 Model 3 - SDE

A_lin3 = A_lin2;
B_lin3 = B_lin2;
E_lin3 = E_lin2;
C_lin3 = C_lin2;
D_lin3 = D_lin2;

% Brownian disturbance intensities from Problem 2.3
sigma3 = 0.5;
sigma4 = 0.5;

Gc = E_lin2*diag([sigma3 sigma4]);


sysd1 = ss(A_lin1,E_lin1,C_lin1,zeros(2,2));
Gdist1 = minreal(tf(sysd1));

disp('Model 1: disturbance -> output')
Gdist1

%% Model 3 transfer functions

% Manipulated inputs -> outputs
sysc3_u = ss(A_lin3,B_lin3,C_lin3,D_lin3);
Glin3 = minreal(tf(sysc3_u));

% Brownian disturbances -> outputs
sysc3_w = ss(A_lin3,Gc,C_lin3,zeros(2,2));
Hnoise3 = minreal(tf(sysc3_w));

disp('Model 3: input -> output')
Glin3

disp('Model 3: Brownian noise -> output')
Hnoise3


%% Model 2 disturbance transfer functions

sysd2 = ss(A_lin2,E_lin2,C_lin2,zeros(2,2));
Gdist2 = minreal(tf(sysd2));

disp('Model 2: disturbance -> output')
Gdist2



%% 5.4 Discrete-time models

Ts = 10;   % sampling time [s]

%% Model 1

% Combine manipulated inputs and disturbances
Bc1 = [B_lin1 E_lin1];

% Continuous-time system
sysc1_all = ss(A_lin1,Bc1,C_lin1,zeros(2,4));

% Zero-order-hold discretization
sysd1 = c2d(sysc1_all,Ts,'zoh');

% Extract discrete-time matrices
[Ad1,Ball_d1,Cd1,Dall_d1] = ssdata(sysd1);

% Separate manipulated inputs and disturbances again
Bd1 = Ball_d1(:,1:2);
Ed1 = Ball_d1(:,3:4);

Ad1
Bd1
Ed1
Cd1

%% Model 2

Bc2 = [B_lin2 E_lin2];

sysc2_all = ss(A_lin2,Bc2,C_lin2,zeros(2,4));

sysd2 = c2d(sysc2_all,Ts,'zoh');

[Ad2,Ball_d2,Cd2,Dall_d2] = ssdata(sysd2);

Bd2 = Ball_d2(:,1:2);
Ed2 = Ball_d2(:,3:4);

Ad2
Bd2
Ed2
Cd2


%% Model 3 - SDE

sysc3_det = ss(A_lin3,B_lin3,C_lin3,D_lin3);

sysd3 = c2d(sysc3_det,Ts,'zoh');

[Ad3,Bd3,Cd3,Dd3] = ssdata(sysd3);


% Continuous diffusion covariance
Qc = Gc*Gc';

n = size(A_lin3,1);

% Van Loan matrix
M = [-A_lin3, Qc;
      zeros(n), A_lin3'] * Ts;

Phi = expm(M);

Phi12 = Phi(1:n,n+1:2*n);
Phi22 = Phi(n+1:2*n,n+1:2*n);

% discrete-time state matrix
Ad3_check = Phi22';

% discrete process-noise covariance
Qd3 = Ad3_check*Phi12;

Ad3
Bd3
Cd3
Qd3


%% 5.5 Markov parameters

Nm = 60;

% Model 1
H1 = zeros(2,2,Nm+1);
H1(:,:,1) = zeros(2,2);     % H0 = D

for k = 1:Nm
    H1(:,:,k+1) = Cd1*(Ad1^(k-1))*Bd1;
end

% Model 2
H2 = zeros(2,2,Nm+1);
H2(:,:,1) = zeros(2,2);

for k = 1:Nm
    H2(:,:,k+1) = Cd2*(Ad2^(k-1))*Bd2;
end

% Model 3
H3 = zeros(2,2,Nm+1);
H3(:,:,1) = zeros(2,2);

for k = 1:Nm
    H3(:,:,k+1) = Cd3*(Ad3^(k-1))*Bd3;
end

tMarkov = (0:Nm)*Ts;


figure

subplot(2,2,1)
stem(tMarkov,squeeze(H1(1,1,:)),'filled')
title('F_1 to h_1')
xlabel('Time [s]')
ylabel('H_{11}[k]')
grid on

subplot(2,2,2)
stem(tMarkov,squeeze(H1(1,2,:)),'filled')
title('F_2 to h_1')
xlabel('Time [s]')
ylabel('H_{12}[k]')
grid on

subplot(2,2,3)
stem(tMarkov,squeeze(H1(2,1,:)),'filled')
title('F_1 to h_2')
xlabel('Time [s]')
ylabel('H_{21}[k]')
grid on

subplot(2,2,4)
stem(tMarkov,squeeze(H1(2,2,:)),'filled')
title('F_2 to h_2')
xlabel('Time [s]')
ylabel('H_{22}[k]')
grid on

sgtitle('Markov parameters - linearized Model 1')

%% 5.5 Compare Markov parameters with Problem 4

P4 = load('problem4_markov.mat');

H_P4 = P4.H;
t_P4 = P4.tMarkov;

%% 5.5 Compare Model 1 with Problem 4

P4 = load('problem4_markov.mat');

H_P4 = P4.H;
t_P4 = P4.tMarkov;

figure

subplot(2,2,1)
plot(t_P4,squeeze(H_P4(1,1,:)),'o-')
hold on
plot(tMarkov,squeeze(H1(1,1,:)),'LineWidth',1.5)
title('F_1 to h_1')
xlabel('Time [s]')
ylabel('H_{11}[k]')
legend('Problem 4 identified','Problem 5 linearized')
grid on

subplot(2,2,2)
plot(t_P4,squeeze(H_P4(1,2,:)),'o-')
hold on
plot(tMarkov,squeeze(H1(1,2,:)),'LineWidth',1.5)
title('F_2 to h_1')
xlabel('Time [s]')
ylabel('H_{12}[k]')
legend('Problem 4 identified','Problem 5 linearized')
grid on

subplot(2,2,3)
plot(t_P4,squeeze(H_P4(2,1,:)),'o-')
hold on
plot(tMarkov,squeeze(H1(2,1,:)),'LineWidth',1.5)
title('F_1 to h_2')
xlabel('Time [s]')
ylabel('H_{21}[k]')
legend('Problem 4 identified','Problem 5 linearized')
grid on

subplot(2,2,4)
plot(t_P4,squeeze(H_P4(2,2,:)),'o-')
hold on
plot(tMarkov,squeeze(H1(2,2,:)),'LineWidth',1.5)
title('F_2 to h_2')
xlabel('Time [s]')
ylabel('H_{22}[k]')
legend('Problem 4 identified','Problem 5 linearized')
grid on

sgtitle('Markov parameters: identified vs linearized Model 1')

Dd1 = zeros(2,2);
Dd2 = zeros(2,2);

save('problem5_statespace.mat', ...
    'Ad1','Bd1','Cd1','Dd1','Ed1', ...
    'Ad2','Bd2','Cd2','Dd2','Ed2', ...
    'Ad3','Bd3','Cd3','Dd3','Qd3','Ts')