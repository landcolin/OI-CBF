function [phi,aux] = backupFlow(p,x,kind,up)
%BACKUPFLOW Integrate paper (8), (12), (29).
if nargin<3, kind='state'; end
if nargin<4, up=oicbf.control(p,x,'primary'); end
x=x(:); n=p.n; t=p.backupGrid;
[f,g]=oicbf.dynamics(p,x); ub=oicbf.control(p,x,'backup');
qp0=f+g*up; qb0=f+g*ub;
if strcmp(p.model,'double_integrator')
    % Exact solution of the same ODE; avoids solver-dependent roundoff.
    phi=[x(1)+x(2)*t-t.^2/2; x(2)-t];
    S=repmat(eye(2),1,1,numel(t)); S(1,2,:)=reshape(t,1,1,[]);
    if strcmp(kind,'stm')
        aux=S;
    elseif strcmp(kind,'directional')
        aux.qp=[qp0(1)+t*qp0(2); repmat(qp0(2),1,numel(t))];
        aux.qb=[qb0(1)+t*qb0(2); repmat(qb0(2),1,numel(t))];
    else
        aux=[];
    end
    return
end
switch kind
    case 'stm', S0=eye(n); z0=[x;S0(:)];
    case 'directional', z0=[x;qp0;qb0];
    case 'state', z0=x;
    otherwise, error('oicbf:BadFlowKind','Unknown backup flow kind.');
end
options=odeset('RelTol',p.rtol,'AbsTol',p.atol);
[~,Y]=ode45(@rhs,t,z0,options); Y=Y';
assert(size(Y,2)==numel(t) && all(isfinite(Y(:))), ...
    'oicbf:IntegrationFailure','Backup integration did not finish with finite values.');
phi=Y(1:n,:);
if strcmp(kind,'stm')
    aux=reshape(Y(n+1:end,:),n,n,[]);
elseif strcmp(kind,'directional')
    aux.qp=Y(n+1:2*n,:); aux.qb=Y(2*n+1:end,:);
else
    aux=[];
end
    function dz=rhs(~,z)
        xx=z(1:n); [ff,gg]=oicbf.dynamics(p,xx);
        dx=ff+gg*oicbf.control(p,xx,'backup');
        if strcmp(kind,'state'), dz=dx; return; end
        F=oicbf.jacobian(p,xx);
        if strcmp(kind,'stm')
            dS=F*reshape(z(n+1:end),n,n); dz=[dx;dS(:)];
        else
            dz=[dx;F*z(n+1:2*n);F*z(2*n+1:end)];
        end
    end
end
