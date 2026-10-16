function [Gm, Gh] = fourTankTransferFunctions(model)
%FOURTANKTRANSFERFUNCTIONS Lower-tank mass and level responses to pump flows.
Gm = minreal(tf(ss(model.A, model.B, model.Cmass, model.D)), 1e-8);
Gh = minreal(tf(ss(model.A, model.B, model.C, model.D)), 1e-8);
end
