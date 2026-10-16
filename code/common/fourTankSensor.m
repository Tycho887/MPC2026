function y = fourTankSensor(x, p, noise)
%FOURTANKSENSOR Measured lower-tank levels [cm], with supplied additive noise.
% Random sampling belongs to the experiment, not the observation model.
if nargin < 3, noise = 0; end
y = fourTankOutput(x, p) + noise;
end
