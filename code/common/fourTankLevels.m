function levels = fourTankLevels(x, p)
%FOURTANKLEVELS Convert 4-by-N mass states [g] to tank levels [cm].
levels = x ./ (p.rho * p.A);
end
