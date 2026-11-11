function [h,dh,hb,dhb,hi,dhi] = barriers(p,x)
%BARRIERS h, terminal h_b and their row gradients: (5), (31), (33).
x=x(:);
if strcmp(p.model,'double_integrator')
    h=-x(1); dh=[-1 0]; hb=-x(2); dhb=[0 -1];
    hi=hb; dhi=dhb; return
end
h=p.ng'*(x(4:5)-p.pg); dh=zeros(1,8); dh(4:5)=p.ng';
hi=[p.cphi^2-(p.phistar-x(1))^2; p.ctheta^2-x(2)^2;
    1-((p.Hstar-x(6))/p.cH)^2; p.cP^2-x(7)^2;
    p.cNz^2-(p.Nzstar-x(8))^2;
    h+p.rho*(p.ng'*[-sin(x(3));cos(x(3))]-1)-p.c6];
dhi=zeros(6,8);
dhi(1,1)=2*(p.phistar-x(1)); dhi(2,2)=-2*x(2);
dhi(3,6)=2*(p.Hstar-x(6))/p.cH^2;
dhi(4,7)=-2*x(7); dhi(5,8)=2*(p.Nzstar-x(8));
dhi(6,:)=dh; dhi(6,3)=p.rho*p.ng'*[-cos(x(3));-sin(x(3))];
if ~p.altitudeNormalized
    hi(3)=p.cH^2-(p.Hstar-x(6))^2;
    dhi(3,6)=2*(p.Hstar-x(6));
end
if strcmp(p.radiusConvention,'signed-corrected')
    hi(6)=hi(6)+p.rho-abs(p.rho);
end
% Stable log-sum-exp and weights, mathematically the printed soft minimum.
z=-p.kappa*hi; zmax=max(z); w=exp(z-zmax); w=w/sum(w);
hb=-(zmax+log(sum(exp(z-zmax))))/p.kappa;
dhb=w'*dhi;
end
