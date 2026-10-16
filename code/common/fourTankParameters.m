function p = fourTankParameters()
%FOURTANKPARAMETERS Physical parameters of the modified four-tank system.
% Units: cm, seconds, grams; input and disturbance flows are cm^3/s.
p.a = 1.2272 * ones(4,1);     % Outlet areas [cm^2]
p.A = 380.1327 * ones(4,1);   % Tank areas [cm^2]
p.gamma = [0.6; 0.7];        % Pump flow fractions to lower tanks
p.g = 981;                  % Gravity [cm/s^2]
p.rho = 1;                  % Water density [g/cm^3]
end
