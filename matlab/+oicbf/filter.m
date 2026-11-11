function [u,info] = filter(p,x,up,mode)
%FILTER Safety filtering: (13)-(15), (18)-(21), (28)-(29).
mode=validatestring(mode,{'bCBF','blending','OI','nominal'});
up=up(:); info=struct('mu',NaN,'feasible',true,'fallback',false, ...
    'solverTime',0,'backupTime',0,'minResidual',NaN,'status','ok');
if strcmp(mode,'nominal'), u=up; return; end
ticBackup=tic;
if strcmp(mode,'blending')
    [phi,~]=oicbf.backupFlow(p,x,'state',up);
    hs=zeros(1,size(phi,2));
    for j=1:size(phi,2), hs(j)=oicbf.barriers(p,phi(:,j)); end
    delta=0;
    % Python applies tightening in blending only for the aircraft example.
    if p.tighteningEnabled && strcmp(p.model,'aircraft')
        delta=.5*p.backupDt*p.Lh*p.supFcl;
    end
    [~,~,hb]=oicbf.barriers(p,phi(:,end));
    hI=min([hs-delta,hb]); info.mu=exp(-p.beta*max(hI,0));
    info.backupTime=toc(ticBackup);
    u=up+info.mu*(oicbf.control(p,x,'backup')-up); info.phi=phi; return
end
[C,d,phi]=oicbf.constraintSystem(p,x,up,mode);
info.backupTime=toc(ticBackup); info.phi=phi; info.C=C; info.d=d;
ticSolver=tic;
if strcmp(mode,'OI')
    positive=C>0;
    info.mu=max(d(positive)./C(positive)); % Paper (21), same formula as Python.
    residual=C*info.mu-d; info.minResidual=min(residual);
    info.feasible=isfinite(info.mu) && info.mu>=-1e-10 && info.mu<=1+1e-10 ...
        && all(isfinite(residual)) && all(residual>=-1e-7);
    u=up+info.mu*(oicbf.control(p,x,'backup')-up);
else
    [u,info.feasible,info.exitflag]=oicbf.solveQP(up,C,d,p.lb,p.ub,p.qpSolver);
    info.qpSolver=p.qpSolver;
    info.minResidual=min(C*u-d);
end
info.solverTime=toc(ticSolver);
if ~info.feasible
    info.status='infeasible';
    switch p.failurePolicy
        case 'error'
            error('oicbf:Infeasible','%s constraints infeasible; do not claim a safe QP solution.',mode);
        case 'backup'
            u=oicbf.control(p,x,'backup'); info.fallback=true;
        case 'legacy'
            % Reproduce Python's unchecked OI value or bCBF backup fallback.
            if strcmp(mode,'bCBF'), u=oicbf.control(p,x,'backup'); info.fallback=true; end
        otherwise, error('oicbf:BadPolicy','Unknown failure policy.');
    end
end
end
