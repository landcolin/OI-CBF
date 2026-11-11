function report = run_batch_tests(samplesPerModel,execution,workers)
%RUN_BATCH_TESTS Deterministic short-horizon tests of varied initial states.
% Each model/sample is run with bCBF, blending and OI, for 40 held-input steps.
% Includes known counterexamples: infeasibility is a recorded outcome.
if nargin<1, samplesPerModel=4; end
if nargin<2, execution='parallel'; end
if nargin<3, workers=2; end
validateattributes(samplesPerModel,{'numeric'},{'scalar','integer','positive'});
validateattributes(workers,{'numeric'},{'scalar','integer','positive'});
execution=validatestring(execution,{'serial','parallel'});
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'test_results','matlab'); if ~isfolder(out), mkdir(out); end
savedRng=rng; restoreRng=onCleanup(@() rng(savedRng)); rng(20261108,'twister');
cases=struct('model',{},'mode',{},'sample',{},'x0',{});
for models={'double_integrator','aircraft'}
    model=models{1}; p=oicbf.config(model);
    for sample=1:samplesPerModel
        x=p.x0;
        if sample==2
            if p.n==2, x=[-.128;.5];
            else, x=[p.phistar;0;0;24990;0;p.Hstar;0;p.Nzstar]; end
        elseif sample>2
            if p.n==2, x=[-2+rand;.2+.4*rand];
            else, x=x+(2*rand(8,1)-1).*[.1;.05;.1;1000;500;100;.05;.1]; end
        end
        for modes={'bCBF','blending','OI'}
            cases(end+1)=struct('model',model,'mode',modes{1},'sample',sample,'x0',x); %#ok<AGROW>
        end
    end
end
poolStartup=0; pool=[]; ownedPool=false;
if strcmp(execution,'parallel')
    assert(~isempty(which('parpool')),'oicbf:MissingParallelToolbox','Install Parallel Computing Toolbox first.');
    pool=gcp('nocreate');
    if isempty(pool)
        timer=tic; pool=parpool('Processes',workers); poolStartup=toc(timer); ownedPool=true;
    end
end
closeOwnedPool=onCleanup(@() closePool(pool,ownedPool)); results=cell(numel(cases),1); timer=tic;
if strcmp(execution,'parallel')
    parfor k=1:numel(cases), results{k}=oneCase(cases(k)); end
    usedWorkers=pool.NumWorkers;
else
    for k=1:numel(cases), results{k}=oneCase(cases(k)); end
    usedWorkers=0;
end
report=struct('execution',execution,'samplesPerModel',samplesPerModel,'stepsPerCase',40, ...
    'workers',usedWorkers,'poolStartupSeconds',poolStartup,'computeSeconds',toc(timer));
report.results=[results{:}];
report.completed=sum(strcmp({report.results.status},'completed'));
report.infeasible=sum(strcmp({report.results.status},'constraint_infeasible'));
report.otherErrors=sum(strcmp({report.results.status},'error'));
report.safetyViolations=sum([report.results.minSafetyMargin]<-1e-8);
report.inputViolations=sum([report.results.maxInputViolation]>1e-8);
fid=fopen(fullfile(out,['batch_' execution '.json']),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
fprintf('%s: %d cases, %d completed, %d infeasible, %d other errors.\n',execution,numel(cases),report.completed,report.infeasible,report.otherErrors);
end

function closePool(pool,ownedPool)
if ownedPool && ~isempty(pool) && isvalid(pool), delete(pool); end
end

function result=oneCase(c)
p=oicbf.config(c.model); p.x0=c.x0; p.steps=40;
result=struct('model',c.model,'mode',c.mode,'sample',c.sample,'x0',c.x0, ...
    'status','completed','minSafetyMargin',NaN,'maxInputViolation',NaN, ...
    'finite',false,'finalState',[],'errorId','','errorMessage','');
try
    r=oicbf.simulate(p,c.mode); ua=r.u(:,2:end);
    result.minSafetyMargin=min(r.h);
    result.maxInputViolation=max([0;reshape(p.lb-ua,[],1);reshape(ua-p.ub,[],1)]);
    result.finite=all(isfinite([r.x(:);r.u(:)])); result.finalState=r.x(:,end);
catch ex
    result.status='error';
    if strcmp(ex.identifier,'oicbf:Infeasible'), result.status='constraint_infeasible'; end
    result.errorId=ex.identifier; result.errorMessage=ex.message;
end
end
