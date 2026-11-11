function [u,du] = control(p,x,kind)
%CONTROL Primary and backup controllers; paper (13), (16), (32).
if nargin<3, kind='primary'; end
x=x(:); du=zeros(p.m,p.n);
if strcmp(p.model,'double_integrator')
    if strcmp(kind,'backup'), u=-1; else, u=1; end
    return
end
if strcmp(kind,'backup')
    phid=p.phistar; Htarget=p.Hstar;
else
    e=p.pPrimary-x(4:5); los=atan2(e(2),e(1));
    epsi=atan2(sin(los-x(3)),cos(los-x(3)));
    if p.absoluteHeading, epsi=abs(epsi); end
    phid=p.phimax*tanh(p.VT/p.g*p.kpsi*epsi/p.phimax);
    Htarget=p.Hprimary;
end
Nzdes=1/cos(phid);
raw=[x(7)+p.tp*(p.kphi*(phid-x(1))-p.kP*x(7));
     x(8)+p.tz*(p.kN*(Nzdes-x(8))-p.ktheta*x(2)+p.kH*(Htarget-x(6)))];
if strcmp(kind,'backup')
    [u1,s1]=oicbf.softSaturate(raw(1),p.lb(1),p.ub(1),true);
    [u2,s2]=oicbf.softSaturate(raw(2),p.lb(2),p.ub(2),false);
    u=[u1;u2];
    du(1,1)=-s1*p.tp*p.kphi; du(1,7)=s1*(1-p.tp*p.kP);
    du(2,2)=-s2*p.tz*p.ktheta; du(2,6)=-s2*p.tz*p.kH;
    du(2,8)=s2*(1-p.tz*p.kN);
else
    u=min(max(raw,p.lb),p.ub);
end
end
