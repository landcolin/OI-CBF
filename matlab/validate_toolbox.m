function report = validate_toolbox()
%VALIDATE_TOOLBOX Independent Optimization Toolbox checks and full bCBF runs.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'test_results','matlab'); if ~isfolder(out), mkdir(out); end
assert(~isempty(which('quadprog')),'oicbf:MissingOptimizationToolbox','Install Optimization Toolbox first.');
data=load(fullfile(root,'test_results','projection_reference.mat'));
report=struct('quadprogPath',which('quadprog'),'cases',numel(data.cases), ...
    'feasibilityMismatches',0,'maxControlErrorVsPython',0,'maxControlErrorVsProjection',0);
for k=1:numel(data.cases)
    c=data.cases{k}; n=double(c.n);
    [u,ok]=oicbf.solveQP(c.up(:),c.C,c.d(:),-ones(n,1),ones(n,1),'quadprog');
    [v,projectionOK]=oicbf.projectPolytope(c.up(:),c.C,c.d(:),-ones(n,1),ones(n,1));
    report.feasibilityMismatches=report.feasibilityMismatches+(ok~=c.feasible)+(ok~=projectionOK);
    if ok && c.feasible
        report.maxControlErrorVsPython=max(report.maxControlErrorVsPython,max(abs(u-c.u(:))));
        report.maxControlErrorVsProjection=max(report.maxControlErrorVsProjection,max(abs(u-v)));
    end
end
report.runs=struct([]); report.maxOIClosedFormDifference=0;
for models={'double_integrator','aircraft'}
    model=models{1}; p=oicbf.config(model); p.qpSolver='quadprog';
    fixture=load(fullfile(root,'test_results',[model '_reference.mat'])); ref=fixture.ref;
    for k=[1 5 9]
        x=ref.x(:,k); up=oicbf.control(p,x); [~,info]=oicbf.filter(p,x,up,'OI');
        [mu,ok]=oicbf.solveQP(0,info.C,info.d,0,1,'quadprog');
        assert(ok,'oicbf:ToolboxOI','OI reference QP unexpectedly failed.');
        report.maxOIClosedFormDifference=max(report.maxOIClosedFormDifference,abs(mu-info.mu));
    end
    fprintf('Full quadprog bCBF simulation: %s\n',model);
    result=oicbf.simulate(p,'bCBF');
    baseline=load(fullfile(out,[model '_bCBF.mat']));
    item=struct('model',model,'minSafetyMargin',min(result.h), ...
        'allQPsFeasible',all(result.feasible), ...
        'maxStateErrorVsProjection',max(abs(result.x-baseline.result.x),[],'all'), ...
        'maxStateErrorVsPython',max(abs(result.x-ref.trajectory_bCBF.x),[],'all'), ...
        'maxControlErrorVsProjection',max(abs(result.u-baseline.result.u),[],'all'), ...
        'maxInputViolation',max([0;reshape(p.lb-result.u(:,2:end),[],1);reshape(result.u(:,2:end)-p.ub,[],1)]));
    item.passed=item.allQPsFeasible && all(isfinite([result.x(:);result.u(:)])) ...
        && item.minSafetyMargin>=-1e-8 && item.maxInputViolation<=1e-8 ...
        && item.maxStateErrorVsProjection<=1e-5 && item.maxControlErrorVsProjection<=1e-6;
    if isempty(report.runs), report.runs=item; else, report.runs(end+1)=item; end
    save(fullfile(out,[model '_quadprog_bCBF.mat']),'result','p');
end
report.passed=report.feasibilityMismatches==0 && report.maxControlErrorVsPython<=1e-8 ...
    && report.maxControlErrorVsProjection<=1e-8 && report.maxOIClosedFormDifference<=1e-8 ...
    && all([report.runs.passed]);
fid=fopen(fullfile(out,'toolbox_validation.json'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); disp(report);
assert(report.passed,'oicbf:ToolboxValidation','Optimization Toolbox validation failed.');
end
