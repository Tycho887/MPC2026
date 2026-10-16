function op = fourTankSteadyState(u, d, p)
%FOURTANKSTEADYSTATE Analytical equilibrium for constant pump/disturbance flows.
[Q, U, D] = fourTankRouting(p);
op.u = u; op.d = d;
op.q = -Q \ (U*u + D*d);
assert(all(op.q >= 0), 'Equilibrium requires nonnegative tank discharge.');
op.h = (op.q ./ p.a).^2 / (2*p.g);
op.x = p.rho * p.A .* op.h;
op.z = fourTankOutput(op.x, p);
op.timeConstants = 2*p.A .* op.h ./ op.q;
end
