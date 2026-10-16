function model = linearizeFourTank(op, p)
%LINEARIZEFOURTANK Jacobians for deviations about a positive equilibrium.
% Pure numerical matrices: no Control System Toolbox is needed here.
assert(all(op.x > 0), 'Linearization requires positive equilibrium masses.');
[Q, U, D] = fourTankRouting(p);
model.A = p.rho * Q * diag(op.q ./ (2*op.x));
model.B = p.rho * U;
model.E = p.rho * D;
model.Cmass = [eye(2), zeros(2)];
model.C = diag(1 ./ (p.rho*p.A(1:2))) * model.Cmass;
model.D = zeros(2);
model.op = op;
end
