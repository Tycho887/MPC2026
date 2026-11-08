function result = markovParameters(G, Ts, count)
%MARKOVPARAMETERS H_0=D; H_k=C*A^(k-1)*B for k>=1 (no impulse/Ts scaling).
% Continuous models are ZOH-discretized; discrete realizations are preserved.
validateattributes(Ts, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'Ts');
validateattributes(count, {'numeric'}, {'real','finite','scalar','integer','positive'}, mfilename, 'count');
if G.Ts == 0
    sys = minreal(ss(c2d(G,Ts,'zoh')),1e-8);
elseif G.Ts > 0 && abs(G.Ts-Ts) <= 10*eps(max(G.Ts,Ts))
    sys = ss(G);
else
    error('markovParameters:SamplingTimeMismatch', ...
        'The discrete model must have the requested sampling time Ts = %g s.',Ts);
end
H = zeros(size(sys.C,1),size(sys.B,2),count);
H(:,:,1) = sys.D;
power = eye(size(sys.A));
for k = 2:count
    H(:,:,k) = sys.C * power * sys.B;
    power = power * sys.A;
end
result = struct('H',H,'t',(0:count-1)'*Ts,'Ts',Ts,'system',sys);
end
