function controller = tuneFourTankPI(G)
%TUNEFOURTANKPI Static decoupling and diagonal PI tuning in G's output units.
controller.W = dcgain(G) \ eye(2);
effective = minreal(G * controller.W, 1e-8);
options = pidtuneOptions('DesignFocus','reference-tracking');
controller.C1 = pidtune(effective(1,1),'PI',options);
controller.C2 = pidtune(effective(2,2),'PI',options);
controller.K = controller.W * [controller.C1 tf(0); tf(0) controller.C2];
controller.closedLoop = feedback(G*controller.K,eye(2));
end
