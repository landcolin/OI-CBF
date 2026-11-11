function result = simulate(p,mode)
%SIMULATE Zero-order-held controls, matching Python array indexing/time span.
if nargin<2, mode='OI'; end
N=p.steps+1; x=zeros(p.n,N); u=zeros(p.m,N); ud=u;
x(:,1)=p.x0; mu=nan(1,p.steps); feasible=true(1,p.steps);
solverTimes=zeros(1,p.steps); backupTimes=solverTimes; interventions=[];
backupTrajs={}; options=odeset('RelTol',p.rtol,'AbsTol',p.atol);
for k=1:p.steps
    ud(:,k+1)=oicbf.control(p,x(:,k),'primary');
    [u(:,k+1),info]=oicbf.filter(p,x(:,k),ud(:,k+1),mode);
    mu(k)=info.mu; feasible(k)=info.feasible;
    solverTimes(k)=info.solverTime; backupTimes(k)=info.backupTime;
    if norm(u(:,k+1)-ud(:,k+1))>=1e-4, interventions(end+1)=k; end %#ok<AGROW>
    if isfield(info,'phi') && mod(k-1,p.backupSaveN)==0
        backupTrajs{end+1}=info.phi; %#ok<AGROW>
    end
    if strcmp(p.model,'double_integrator')
        x(:,k+1)=[x(1,k)+p.dt*x(2,k)+.5*p.dt^2*u(k+1);x(2,k)+p.dt*u(k+1)];
    else
        [~,xx]=ode45(@(~,xx) heldDynamics(xx,u(:,k+1)),[0 p.dt],x(:,k),options);
        x(:,k+1)=xx(end,:)';
    end
    assert(all(isfinite(x(:,k+1))),'oicbf:NonfiniteState','State integration became nonfinite.');
end
h=zeros(1,N);
for k=1:N, h(k)=oicbf.barriers(p,x(:,k)); end
result=struct('model',p.model,'profile',p.profile,'mode',mode,'t',(0:p.steps)*p.dt, ...
    'x',x,'u',u,'uDesired',ud,'mu',mu,'h',h,'feasible',feasible, ...
    'interventions',interventions,'solverTimes',solverTimes,'backupTimes',backupTimes, ...
    'avgSolverTime',mean(solverTimes),'maxSolverTime',max(solverTimes));
result.backupTrajs=backupTrajs;
    function dx=heldDynamics(xx,uu)
        [f,g]=oicbf.dynamics(p,xx); dx=f+g*uu;
    end
end
