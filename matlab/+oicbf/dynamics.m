function [f,g] = dynamics(p,x)
%DYNAMICS Paper (1), (16), (30); column-vector state convention.
x = x(:);
if strcmp(p.model,'double_integrator')
    f = [x(2);0]; g = [0;1]; return
end
ph=x(1); th=x(2); ps=x(3); P=x(7); Nz=x(8);
f = [P+Nz*p.g/p.VT*sin(ph)*tan(th);
     p.g/p.VT*(Nz*cos(ph)-cos(th));
     Nz*p.g*sin(ph)/(p.VT*cos(th));
     p.VT*cos(th)*cos(ps); p.VT*cos(th)*sin(ps); p.VT*sin(th);
     -P/p.tp; -Nz/p.tz];
g=zeros(8,2); g(7,1)=1/p.tp; g(8,2)=1/p.tz;
end
