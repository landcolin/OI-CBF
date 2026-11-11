function J = jacobian(p,x)
%JACOBIAN Analytic derivative of f_cl=f+g*k_b, paper (12), (29).
if strcmp(p.model,'double_integrator'), J=[0 1;0 0]; return; end
ph=x(1); th=x(2); ps=x(3); Nz=x(8); a=p.g/p.VT;
J=zeros(8);
J(1,1)=Nz*a*cos(ph)*tan(th);
J(1,2)=Nz*a*sin(ph)/cos(th)^2; J(1,7)=1;
J(1,8)=a*sin(ph)*tan(th);
J(2,1)=-Nz*a*sin(ph); J(2,2)=a*sin(th); J(2,8)=a*cos(ph);
J(3,1)=Nz*a*cos(ph)/cos(th);
J(3,2)=Nz*a*sin(ph)*sin(th)/cos(th)^2;
J(3,8)=a*sin(ph)/cos(th);
J(4,2)=-p.VT*sin(th)*cos(ps); J(4,3)=-p.VT*cos(th)*sin(ps);
J(5,2)=-p.VT*sin(th)*sin(ps); J(5,3)=p.VT*cos(th)*cos(ps);
J(6,2)=p.VT*cos(th); J(7,7)=-1/p.tp; J(8,8)=-1/p.tz;
[~,g]=oicbf.dynamics(p,x); [~,du]=oicbf.control(p,x,'backup');
J=J+g*du;
end
