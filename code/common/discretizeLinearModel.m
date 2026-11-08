function discrete = discretizeLinearModel(model, Ts, method)
%DISCRETIZELINEARMODEL Exact ZOH for a continuous numeric model with held disturbances.
% method is 'zoh' (Control System Toolbox) or 'expm' (base MATLAB).
% Copies model metadata and output matrices; replaces A, B, E and sets Ts.
validateattributes(Ts, {'numeric'}, {'real','finite','scalar','positive'}, mfilename, 'Ts');
method = validatestring(method, {'zoh','expm'}, mfilename, 'method');
nx = size(model.A,1);
nu = size(model.B,2);
nd = size(model.E,2);
inputs = [model.B, model.E];

switch method
    case 'zoh'
        feedthrough = [model.D, zeros(size(model.C,1),nd)];
        system = c2d(ss(model.A,inputs,model.C,feedthrough),Ts,'zoh');
        Abar = system.A;
        inputsBar = system.B;
    case 'expm'
        % The augmented exponential also works when A is singular.
        Phi = expm(Ts * [model.A, inputs; zeros(nu+nd,nx+nu+nd)]);
        Abar = Phi(1:nx,1:nx);
        inputsBar = Phi(1:nx,nx+1:nx+nu+nd);
end

discrete = model;
discrete.A = Abar;
discrete.B = inputsBar(:,1:nu);
discrete.E = inputsBar(:,nu+1:nu+nd);
discrete.Ts = Ts;
end
