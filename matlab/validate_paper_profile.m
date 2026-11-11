function report = validate_paper_profile()
%VALIDATE_PAPER_PROFILE Compare implementation with independently written
% literal paper barriers and STM chain-rule construction of (19)-(20).
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'test_results','matlab'); if ~isfolder(out), mkdir(out); end
report=struct; report.models=struct([]);
for names={'double_integrator','aircraft'}
    name=names{1}; p=oicbf.config(name,'paper'); p.rtol=1e-10; p.atol=1e-10;
    data=load(fullfile(root,'test_results',[name '_reference.mat'])); ref=data.ref;
    barrierError=0; gradientError=0; coeffError=0;
    for k=[1 5 9]
        x=ref.x(:,k); u=oicbf.control(p,x,'primary');
        [h,dh,hb,dhb,hi]=oicbf.barriers(p,x);
        [expected,components]=printedBarriers(p,x);
        barrierError=max(barrierError,max(abs([h;hb;hi]-[expected;components])));
        D=numericGradient(@(xx) printedBarriers(p,xx),x);
        gradientError=max(gradientError,max(abs([dh;dhb]-D),[],'all'));
        [C,d]=oicbf.constraintSystem(p,x,u,'OI');
        [phi,S]=oicbf.backupFlow(p,x,'stm',u);
        [f,g]=oicbf.dynamics(p,x); fp=f+g*u; fb=f+g*oicbf.control(p,x,'backup');
        Cr=zeros(size(C)); dr=zeros(size(d)); Cr(1:2)=[1;-1]; dr(1:2)=[0;-1];
        % Independent route: full STM (12), numeric barrier gradients,
        % then chain rule (11), instead of directional integration (29).
        for j=1:size(phi,2)
            val=printedBarriers(p,phi(:,j)); grad=numericGradient(@(xx) printedBarriers(p,xx),phi(:,j));
            hp=grad(1,:)*S(:,:,j)*fp; hback=grad(1,:)*S(:,:,j)*fb;
            Cr(j+2)=hback-hp; dr(j+2)=-hp-p.alpha(val(1));
        end
        val=printedBarriers(p,phi(:,end)); grad=numericGradient(@(xx) printedBarriers(p,xx),phi(:,end));
        hp=grad(2,:)*S(:,:,end)*fp; hback=grad(2,:)*S(:,:,end)*fb;
        Cr(end)=hback-hp; dr(end)=-hp-p.alphaB(val(2));
        coeffError=max(coeffError,max(abs([C-Cr;d-dr])));
    end
    item=struct('model',name,'states',3,'maxBarrierError',barrierError, ...
        'maxGradientError',gradientError,'maxOIConstraintCoefficientError',coeffError, ...
        'barrierTolerance',1e-9,'gradientTolerance',1e-5,'coefficientTolerance',1e-5);
    item.passed=barrierError<=item.barrierTolerance && gradientError<=item.gradientTolerance ...
        && coeffError<=item.coefficientTolerance;
    if isempty(report.models), report.models=item; else, report.models(end+1)=item; end
end
report.allPassed=all([report.models.passed]);
report.interpretation='Literal formula implementation verified numerically; paper safety claims are not established.';
fid=fopen(fullfile(out,'paper_profile_validation.json'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); disp(report.models);
assert(report.allPassed,'oicbf:PaperProfileValidation','Paper-profile formula check failed.');
end

function [values,hi]=printedBarriers(p,x)
if p.n==2, values=[-x(1);-x(2)]; hi=-x(2); return; end
h=p.ng'*(x(4:5)-p.pg);
rho=p.VT^2/(p.g*tan(p.phistar));
hi=[p.cphi^2-(p.phistar-x(1))^2;p.ctheta^2-x(2)^2;
    p.cH^2-(p.Hstar-x(6))^2;p.cP^2-x(7)^2;
    p.cNz^2-(1/cos(p.phistar)-x(8))^2;
    h+rho*(p.ng'*[-sin(x(3));cos(x(3))]-1)-p.c6];
% Stable evaluation of the literal logarithmic soft minimum.
smallest=min(hi); hb=smallest-log(sum(exp(-p.kappa*(hi-smallest))))/p.kappa;
values=[h;hb];
end
function D=numericGradient(fun,x)
D=zeros(2,numel(x));
for i=1:numel(x)
    step=eps^(1/3)*max(1,abs(x(i))); perturb=zeros(size(x)); perturb(i)=step;
    D(:,i)=(fun(x+perturb)-fun(x-perturb))/(2*step);
end
end
