function z = fourTankOutput(x, p)
%FOURTANKOUTPUT Noise-free controlled levels, tanks 1 and 2 [cm].
levels = fourTankLevels(x, p);
z = levels(1:2,:);
end
