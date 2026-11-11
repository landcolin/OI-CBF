function [u,feasible] = projectPolytope(up,C,d,lb,ub)
%PROJECTPOLYTOPE Exact 1D/2D Euclidean projection for the bCBF objective.
% No Optimization Toolbox needed. 2D: clip the input box by halfspaces,
% then project onto its boundary segments. Equivalent to min ||u-up||^2.
up=up(:); lb=lb(:); ub=ub(:); d=d(:);
tol=1e-10; u=nan(size(up)); feasible=false;
if any(~isfinite([C(:);d;up;lb;ub])), return; end
if isscalar(up)
    lo=lb; hi=ub;
    for i=1:numel(d)
        if C(i)>0
            lo=max(lo,d(i)/C(i));
        elseif C(i)<0
            hi=min(hi,d(i)/C(i));
        elseif d(i)>tol
            return
        end
    end
    feasible=lo<=hi+tol;
    if feasible, u=min(max(up,lo),hi); end
    return
end
assert(numel(up)==2,'oicbf:QPDimension','Projection supports the two supplied models.');
if all(C*up>=d-tol) && all(up>=lb-tol) && all(up<=ub+tol)
    u=up; feasible=true; return
end
poly=[lb(1) lb(2);ub(1) lb(2);ub(1) ub(2);lb(1) ub(2)];
for i=1:numel(d)
    if isempty(poly), return; end
    next=zeros(0,2); a=C(i,:); rhs=d(i);
    for j=1:size(poly,1)
        A=poly(j,:); B=poly(mod(j,size(poly,1))+1,:);
        da=a*A'-rhs; db=a*B'-rhs;
        ina=da>=-tol; inb=db>=-tol;
        if ina, next(end+1,:)=A; end %#ok<AGROW>
        if ina~=inb
            next(end+1,:)=A+(da/(da-db))*(B-A); %#ok<AGROW>
        end
    end
    poly=next;
end
if isempty(poly), return; end
best=inf;
for j=1:size(poly,1)
    a=poly(j,:)'; b=poly(mod(j,size(poly,1))+1,:)'; v=b-a;
    if dot(v,v)>0, t=min(max(dot(up-a,v)/dot(v,v),0),1); else, t=0; end
    candidate=a+t*v; cost=sum((candidate-up).^2);
    if cost<best, best=cost; u=candidate; end
end
feasible=all(C*u>=d-1e-7);
end
