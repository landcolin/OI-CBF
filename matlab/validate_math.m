function report = validate_math(outputDir)
%VALIDATE_MATH Independent formula and Python-reference checks, both models.
% This tests implementation consistency; paper counterexamples are recorded
% explicitly and do not become false claims of a proven safety theorem.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1, outputDir=fullfile(root,'test_results','matlab'); end
if ~isfolder(outputDir), mkdir(outputDir); end
report=struct; report.checks=struct('name',{},'passed',{},'error',{},'tolerance',{});
report.matlabVersion=version;
for models={'double_integrator','aircraft'}
    model=models{1}; p=oicbf.config(model,'source');
    fixture=load(fullfile(root,'test_results',[model '_reference.mat'])); ref=fixture.ref;
    % NumPy's constant A/B and gradients can be integer arrays in MAT files.
    for field={'g','J','dh','dhb','dhi'}, ref.(field{1})=double(ref.(field{1})); end
    errF=0; errG=0; errPrimary=0; errBackup=0; errJ=0; errH=0; errHi=0; errHb=0;
    errGrad=0; errJacFD=0; errGradFD=0; errPaperF=0;
    for k=1:size(ref.x,2)
        x=ref.x(:,k); [f,g]=oicbf.dynamics(p,x);
        [h,dh,hb,dhb,hi,dhi]=oicbf.barriers(p,x);
        errF=max(errF,maxabs(f-ref.f(:,k))); errG=max(errG,maxabs(g-ref.g(:,:,k)));
        errPrimary=max(errPrimary,maxabs(oicbf.control(p,x,'primary')-ref.primary(:,k)));
        errBackup=max(errBackup,maxabs(oicbf.control(p,x,'backup')-ref.backup(:,k)));
        errJ=max(errJ,maxabs(oicbf.jacobian(p,x)-ref.J(:,:,k)));
        errH=max(errH,abs(h-ref.h(k))); errHi=max(errHi,maxabs(hi-ref.hi(:,k)));
        % Independent stable expression based on exported component barriers.
        if p.n==8
            zz=-p.kappa*ref.hi(:,k); z0=max(zz); ww=exp(zz-z0); ww=ww/sum(ww);
            hbRef=-(z0+log(sum(exp(zz-z0))))/p.kappa;
            dhbRef=ww'*ref.dhi(:,:,k);
        else
            hbRef=ref.hb(k); dhbRef=ref.dhb(k,:);
        end
        errHb=max(errHb,abs(hb-hbRef));
        errGrad=max(errGrad,maxabs([dh-ref.dh(k,:);dhb-dhbRef]));
        errPaperF=max(errPaperF,maxabs(f-paperF(p,x)));
        Jfd=finiteDiff(@(xx) backupDynamics(p,xx),x);
        errJacFD=max(errJacFD,maxabs(Jfd-oicbf.jacobian(p,x)));
        gradFD=finiteDiff(@(xx) allBarriers(p,xx),x);
        errGradFD=max(errGradFD,maxabs(gradFD-[dh;dhb;dhi]));
    end
    add([model '/dynamics-python'],max(errF,errG),1e-10);
    add([model '/dynamics-paper-16-or-30'],errPaperF,1e-12);
    add([model '/primary-python'],errPrimary,1e-10);
    add([model '/backup-python'],errBackup,1e-10);
    add([model '/jacobian-python'],errJ,1e-10);
    add([model '/barriers-python-components'],max(errH,errHi),1e-9);
    add([model '/stable-softmin-formula'],errHb,1e-9);
    add([model '/barrier-gradients-reference'],errGrad,1e-9);
    add([model '/jacobian-finite-differences'],errJacFD,2e-6);
    add([model '/barrier-gradients-finite-differences'],errGradFD,5e-5);

    flowState=0; flowS=0; flowQ=0; directionIdentity=0; cErr=0; dErr=0; uErr=0; qpErr=0;
    for k=1:numel(ref.flow)
        ff=ref.flow{k}; x=ff.x(:); up=oicbf.control(p,x,'primary');
        [phi,S]=oicbf.backupFlow(p,x,'stm',up);
        [~,q]=oicbf.backupFlow(p,x,'directional',up);
        flowState=max(flowState,maxabs(phi-ff.phi)); flowS=max(flowS,maxabs(S-ff.S));
        flowQ=max(flowQ,max(maxabs(q.qp-ff.qp),maxabs(q.qb-ff.qb)));
        tight=p; tight.rtol=1e-10; tight.atol=1e-10;
        [~,St]=oicbf.backupFlow(tight,x,'stm',up); [~,qt]=oicbf.backupFlow(tight,x,'directional',up);
        [f,g]=oicbf.dynamics(tight,x); fp=f+g*up; fb=f+g*oicbf.control(tight,x,'backup');
        for j=1:size(St,3)
            directionIdentity=max(directionIdentity,max(maxabs(qt.qp(:,j)-St(:,:,j)*fp),maxabs(qt.qb(:,j)-St(:,:,j)*fb)));
        end
        [C,d]=oicbf.constraintSystem(p,x,up,'OI');
        cErr=max(cErr,maxabs(C-ref.OI_C(:,k))); dErr=max(dErr,maxabs(d-ref.OI_d(:,k)));
        [u,info]=oicbf.filter(p,x,up,'OI'); uErr=max(uErr,maxabs(u-ref.OI_u(:,k)));
        [muQP,ok]=oicbf.projectPolytope(0,C,d,0,1);
        if ~ok, error('oicbf:ReferenceQP','Reference OI constraints unexpectedly infeasible.'); end
        add(sprintf('%s/OI-KKT-%d',model,k),abs(info.mu-muQP),1e-8);
        [uQP,ok]=oicbf.projectPolytope(up,ref.QP_C(:,:,k),ref.QP_d(:,k),p.lb,p.ub);
        if ~ok, error('oicbf:ReferenceQP','Python QP fixture unexpectedly infeasible.'); end
        qpErr=max(qpErr,maxabs(uQP-ref.QP_u(:,k)));
    end
    if p.n==2, flowTolerance=1e-9; else, flowTolerance=.1; end
    add([model '/backup-state-python-default-tolerances'],flowState,flowTolerance);
    add([model '/STM-python-default-tolerances'],flowS,flowTolerance);
    add([model '/directional-python-default-tolerances'],flowQ,flowTolerance);
    add([model '/eq29-vs-eq12-tight-integration'],directionIdentity,1e-4);
    add([model '/OI-coefficients-python-default-tolerances'],max(cErr,dErr),flowTolerance);
    add([model '/OI-control-python-default-tolerances'],uErr,2e-4);
    add([model '/polytope-projection-vs-python-quadprog'],qpErr,1e-7);
end

% Independent tests of the QP solver, including negative/zero coefficients.
[u,ok]=oicbf.projectPolytope([2;2],[1 1],1,[-1;-1],[1;1]);
add('projection/box-corner',maxabs(u-[1;1])+double(~ok),1e-12);
[u,ok]=oicbf.projectPolytope([1;1],[-1 -1],-1,[-2;-2],[2;2]);
add('projection/oblique-edge',maxabs(u-[.5;.5])+double(~ok),1e-12);
[~,ok]=oicbf.projectPolytope(0,[1;-1],[2;-1],0,1);
add('projection/infeasible-interval-rejected',double(ok),0);
[~,ok]=oicbf.projectPolytope(0,0,1,0,1);
add('projection/infeasible-zero-row-rejected',double(ok),0);

p=oicbf.config('double_integrator'); caught=false;
try
    oicbf.filter(p,[-.128;.5],1,'OI');
catch ex
    caught=strcmp(ex.identifier,'oicbf:Infeasible');
end
add('OI/infeasibility-guard',double(~caught),0);
p.failurePolicy='legacy'; [u,info]=oicbf.filter(p,[-.128;.5],1,'OI');
add('OI/explicit-legacy-reproduces-original-bug',abs(u+1.0409090410181818),1e-10);
report.originalBoundaryCounterexample=struct('mu',info.mu,'u',u,'feasible',info.feasible);

% Printed formula differences and the signed-radius counterexample.
source=oicbf.config('aircraft','source'); paper=oicbf.config('aircraft','paper');
x=paper.x0; x(6)=10150;
[~,~,~,~,hSource]=oicbf.barriers(source,x); [~,~,~,~,hPaper]=oicbf.barriers(paper,x);
add('paper/altitude-h3-unnormalized',abs(hPaper(3)-(paper.cH^2-150^2)),1e-12);
add('source/altitude-h3-normalized',abs(hSource(3)-(1-(150/400)^2)),1e-12);
add('paper/terminal-gain-alphaB',abs(paper.alphaOIBackup(.2)-paper.alphaB(.2)),1e-12);
add('source/terminal-gain-alpha',abs(source.alphaOIBackup(.2)-source.alpha(.2)),1e-12);
x=paper.x0; x(3)=.5; uS=oicbf.control(source,x); uP=oicbf.control(paper,x);
add('paper/signed-heading-no-abs',double(~(uS(1)>0 && uP(1)<0)),0);
x=[paper.phistar;0;0;26000;0;paper.Hstar;0;paper.Nzstar];
[h,~,~,~,hi]=oicbf.barriers(paper,x);
report.printedRadiusCounterexample=struct('rho',paper.rho,'h',h,'h6',hi(6), ...
    'satisfiesEq33Equalities',true,'invalidSubsetClaim',h<0 && hi(6)>0);
add('audit/printed-radius-counterexample-reproduced',double(~report.printedRadiusCounterexample.invalidSubsetClaim),0);
paper.radiusConvention='signed-corrected'; [~,~,~,~,hc]=oicbf.barriers(paper,x);
add('amendment/signed-radius-excludes-counterexample',double(hc(6)>=0),0);
for j=[1 3 6 7]
    xx=source.x0;
    if j==1, xx(7)=-30; elseif j==3, xx(6)=12000; elseif j==6, xx(7)=1e8; else, xx(6)=1e8; end
    uu=oicbf.control(source,xx,'backup'); [~,~,hh,dd]=oicbf.barriers(source,xx);
    add(sprintf('stability/finite-%d',j),double(~all(isfinite([uu;hh;dd(:)]))),0);
end
report.allImplementationChecksPassed=all([report.checks.passed]);
report.paperSafetyClaimsValidated=false;
writeJSON(fullfile(outputDir,'math_validation.json'),report);
save(fullfile(outputDir,'math_validation.mat'),'report');
disp(struct2table(report.checks));
assert(report.allImplementationChecksPassed,'oicbf:MathValidation','Some implementation checks failed; see math_validation.json.');
    function add(name,err,tolerance)
        report.checks(end+1)=struct('name',name,'passed',isfinite(err)&&err<=tolerance,'error',err,'tolerance',tolerance);
    end
end

function y=paperF(p,x)
if p.n==2, y=[x(2);0]; return; end
% Literal Eq. (30), kept independent of oicbf.dynamics.
y=[x(7)+x(8)*p.g/p.VT*sin(x(1))*tan(x(2));
   p.g/p.VT*(x(8)*cos(x(1))-cos(x(2)));
   x(8)*p.g*sin(x(1))/(p.VT*cos(x(2)));
   p.VT*cos(x(2))*cos(x(3)); p.VT*cos(x(2))*sin(x(3));
   p.VT*sin(x(2)); -x(7)/p.tp; -x(8)/p.tz];
end
function y=backupDynamics(p,x)
[f,g]=oicbf.dynamics(p,x); y=f+g*oicbf.control(p,x,'backup');
end
function y=allBarriers(p,x)
[h,~,hb,~,hi]=oicbf.barriers(p,x); y=[h;hb;hi];
end
function J=finiteDiff(fun,x)
y=fun(x); J=zeros(numel(y),numel(x));
for k=1:numel(x)
    step=eps^(1/3)*max(1,abs(x(k))); delta=zeros(size(x)); delta(k)=step;
    J(:,k)=(fun(x+delta)-fun(x-delta))/(2*step);
end
end
function e=maxabs(x)
if any(~isfinite(x(:))), e=Inf; else, e=max(abs(x(:))); end
end
function writeJSON(file,value)
fid=fopen(file,'w','n','UTF-8'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
