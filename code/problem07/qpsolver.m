function [x,info] = qpsolver(H,g,l,u,A,bl,bu,xinit)

%constraint bounds as column vectors
bl = bl(:);
bu = bu(:);

%finite upper and lower constraint bounds
iUpper = isfinite(bu);
iLower = isfinite(bl);

%Convert constraints to quadprog format
Aineq = [A(iUpper,:); -A(iLower,:)];
bineq = [bu(iUpper); -bl(iLower)];

%run without displaying output
options = optimoptions('quadprog','Display','off');

[x,~,info] = quadprog(H,g,Aineq,bineq,[],[],l,u,xinit,options);

end
