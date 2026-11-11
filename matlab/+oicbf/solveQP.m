function [u,feasible,exitflag] = solveQP(up,C,d,lb,ub,solver)
%SOLVEQP min 1/2*||u-up||^2 subject to C*u>=d and input bounds.
% MATLAB quadprog uses H=I, f=-up, A=-C, b=-d.
if nargin<6, solver='projection'; end
solver=validatestring(solver,{'projection','quadprog'});
if strcmp(solver,'projection')
    [u,feasible]=oicbf.projectPolytope(up,C,d,lb,ub);
    if feasible, exitflag=1; else, exitflag=-2; end
    return
end
if isempty(which('quadprog'))
    error('oicbf:MissingOptimizationToolbox','quadprog is unavailable. Install Optimization Toolbox or use projection.');
end
persistent options
if isempty(options)
    options=optimoptions('quadprog','Algorithm','active-set','Display','off', ...
        'ConstraintTolerance',1e-10,'OptimalityTolerance',1e-10);
end
up=up(:); lb=lb(:); ub=ub(:); d=d(:);
x0=min(max(up,lb),ub);
[u,~,exitflag]=quadprog(eye(numel(up)),-up,-C,-d,[],[],lb,ub,x0,options);
feasible=exitflag>0 && numel(u)==numel(up) && all(isfinite(u)) ...
    && all(C*u>=d-1e-7) && all(u>=lb-1e-7) && all(u<=ub+1e-7);
if isempty(u), u=nan(size(up)); end
end
