function [C,d,phi,aux] = constraintSystem(p,x,up,mode)
%CONSTRAINTSYSTEM C*u>=d (bCBF) or C*mu>=d (OI): (10),(17),(19),(20).
n=p.n; [f,g]=oicbf.dynamics(p,x);
if strcmp(mode,'OI')
    [phi,aux]=oicbf.backupFlow(p,x,'directional',up); C=[1;-1]; d=[0;-1];
elseif strcmp(mode,'bCBF')
    [phi,aux]=oicbf.backupFlow(p,x,'stm',up);
    C=[eye(p.m);-eye(p.m)]; d=[p.lb;-p.ub];
else
    error('oicbf:BadMode','Constraint mode must be OI or bCBF.');
end
delta=0;
if p.tighteningEnabled, delta=.5*p.backupDt*p.Lh*p.supFcl; end
first=2; if p.includeInitialConstraint, first=1; end
N=size(phi,2); rows=N-first+2; start=size(C,1);
if strcmp(mode,'OI'), C=[C;zeros(rows,1)];
else, C=[C;zeros(rows,p.m)]; end
d=[d;zeros(rows,1)];
for j=first:N
    [h,dh]=oicbf.barriers(p,phi(:,j)); r=start+j-first+1;
    if strcmp(mode,'OI')
        hp=dh*aux.qp(:,j); hb=dh*aux.qb(:,j);
        C(r)=hb-hp; d(r)=-hp-p.alpha(h-delta);
    else
        S=aux(:,:,j); C(r,:)=dh*S*g; d(r)=-dh*S*f-p.alpha(h-delta);
    end
end
[~,~,hb,dhb]=oicbf.barriers(p,phi(:,end));
if strcmp(mode,'OI')
    hp=dhb*aux.qp(:,end); hback=dhb*aux.qb(:,end);
    C(end)=hback-hp; d(end)=-hp-p.alphaOIBackup(hb);
else
    S=aux(:,:,end); C(end,:)=dhb*S*g; d(end)=-dhb*S*f-p.alphaB(hb);
end
assert(numel(x)==n,'oicbf:BadState','Incorrect state dimension.');
end
