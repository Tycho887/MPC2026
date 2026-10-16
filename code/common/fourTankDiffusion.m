function sigma = fourTankDiffusion(p, flowDiffusion)
%FOURTANKDIFFUSION Constant 4-by-2 SDE diffusion matrix [g/s^(1/2)].
% flowDiffusion contains the two disturbance intensities [cm^3/s^(1/2)].
[~, ~, D] = fourTankRouting(p);
sigma = p.rho * D * diag(flowDiffusion);
end
