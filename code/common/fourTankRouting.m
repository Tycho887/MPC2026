function [Q, U, D] = fourTankRouting(p)
%FOURTANKROUTING Drainage, pump, and disturbance routing matrices.
Q = [-1 0 1 0; 0 -1 0 1; 0 0 -1 0; 0 0 0 -1];
U = [p.gamma(1) 0; 0 p.gamma(2); 0 1-p.gamma(2); 1-p.gamma(1) 0];
D = [0 0; 0 0; 1 0; 0 1];
end
