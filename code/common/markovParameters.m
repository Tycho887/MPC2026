function result = markovParameters(G, Ts, count)
%MARKOVPARAMETERS H_0=D; H_k=C*A^(k-1)*B for k>=1 (no impulse/Ts scaling).
sys = minreal(ss(c2d(G,Ts,'zoh')),1e-8);
H = zeros(size(sys.C,1),size(sys.B,2),count);
H(:,:,1) = sys.D;
power = eye(size(sys.A));
for k = 2:count
    H(:,:,k) = sys.C * power * sys.B;
    power = power * sys.A;
end
result = struct('H',H,'t',(0:count-1)'*Ts,'Ts',Ts,'system',sys);
end
