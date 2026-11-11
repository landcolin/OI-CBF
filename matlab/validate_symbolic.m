function report = validate_symbolic()
%VALIDATE_SYMBOLIC Independently differentiate dynamics and component barriers.
% Symbolic expressions use exact decimal constants; simulations stay numeric.
assert(~isempty(which('sym')),'oicbf:MissingSymbolicToolbox','Install Symbolic Math Toolbox first.');
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'test_results','matlab'); if ~isfolder(out), mkdir(out); end
x=sym('x',[8 1],'real'); ph=x(1); th=x(2); ps=x(3); P=x(7); Nz=x(8);
g=sym('9.81'); VT=sym(200); phistar=-sym(pi)/4;
kphi=sym(3); kP=sym(2); kN=sym('.4'); ktheta=sym('.65'); kH=sym('.0001');
raw=[P+kphi*(phistar-ph)-kP*P;
    Nz+kN*(1/cos(phistar)-Nz)-ktheta*th+kH*(10000-x(6))];
A=sym(pi)/2; s=sym(30); epsilon=sym('0.000001');
gamma1=1/(tanh(s*A/2)+epsilon); L=sym(-1); U=sym(4);
offset=(log(1+exp(s*L))-log(1+exp(-s*U)))/s;
slope=1/(1+exp(s*(offset+L)))-1/(1+exp(s*(offset+U)));
gamma2=1/(slope+epsilon);
ub=[-A+(log(1+exp(s*(gamma1*raw(1)+A)))-log(1+exp(s*(gamma1*raw(1)-A))))/s;
    L+(log(1+exp(s*(gamma2*raw(2)-offset-L)))-log(1+exp(s*(gamma2*raw(2)-offset-U))))/s];
% Literal Eq. (30) combined with the source's smooth Eq. (32) controller.
fcl=[P+Nz*g/VT*sin(ph)*tan(th);g/VT*(Nz*cos(ph)-cos(th));
    Nz*g*sin(ph)/(VT*cos(th));VT*cos(th)*cos(ps);VT*cos(th)*sin(ps);
    VT*sin(th);ub(1)-P;ub(2)-Nz];
Js=jacobian(fcl,x);
fprintf('Generating symbolic aircraft Jacobian...\n');
Jfun=matlabFunction(Js,'Vars',{x},'File',fullfile(out,'symbolic_aircraft_jacobian.m'));
previousPath=path; restorePath=onCleanup(@() path(previousPath)); addpath(out); rehash;
p=oicbf.config('aircraft'); data=load(fullfile(root,'test_results','aircraft_reference.mat'));
maxJ=0;
for k=1:size(data.ref.x,2)
    xx=data.ref.x(:,k); delta=Jfun(xx)-oicbf.jacobian(p,xx);
    assert(all(isfinite(delta(:))),'oicbf:SymbolicNonfinite','Generated Jacobian became nonfinite.');
    maxJ=max(maxJ,max(abs(delta(:))));
end
h=25000-x(4); rho=VT^2/(g*tan(phistar));
hi=[(sym(pi)/12)^2-(phistar-ph)^2;(sym(pi)/12)^2-th^2;
    1-((10000-x(6))/400)^2;(sym(pi)/12)^2-P^2;
    sym('.1')^2-(1/cos(phistar)-Nz)^2;h+rho*(sin(ps)-1)-5];
componentFun=matlabFunction(jacobian(hi,x),'Vars',{x});
paperHi=hi; paperHi(3)=sym(400)^2-(10000-x(6))^2;
paperComponentFun=matlabFunction(jacobian(paperHi,x),'Vars',{x});
maxGradient=0;
for profiles={'source','paper'}
    p=oicbf.config('aircraft',profiles{1});
    for k=1:size(data.ref.x,2)
        xx=data.ref.x(:,k); [~,~,~,~,~,actual]=oicbf.barriers(p,xx);
        if strcmp(profiles{1},'source'), expected=componentFun(xx); else, expected=paperComponentFun(xx); end
        maxGradient=max(maxGradient,max(abs(actual-expected),[],'all'));
    end
end
% Exact identity for the soft-min gradient, independent of its stable evaluator.
y=sym('y',[6 1],'real'); kappa=sym('kappa','positive');
weights=exp(-kappa*y)/sum(exp(-kappa*y));
soft=-log(sum(exp(-kappa*y)))/kappa;
softIdentity=isAlways(simplify(jacobian(soft,y)-weights')==0);
z=sym('z',[2 1],'real'); Jdi=jacobian([z(2);-sym(1)],z);
doubleIdentity=isequal(Jdi,sym([0 1;0 0]));
report=struct('matlabVersion',version,'aircraftStates',size(data.ref.x,2), ...
    'maxJacobianError',maxJ,'maxComponentGradientError',maxGradient, ...
    'softminGradientIdentity',all(softIdentity),'doubleIntegratorJacobianIdentity',doubleIdentity);
report.passed=maxJ<1e-10 && maxGradient<1e-9 && report.softminGradientIdentity && doubleIdentity;
report.interpretation='Derivative implementation checks, not a proof of the paper safety assumptions.';
fid=fopen(fullfile(out,'symbolic_validation.json'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); disp(report);
assert(report.passed,'oicbf:SymbolicValidation','Symbolic derivative checks failed.');
end
