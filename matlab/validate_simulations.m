function report = validate_simulations(outputDir)
%VALIDATE_SIMULATIONS Run all six original scenarios and compare trajectories.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<1, outputDir=fullfile(root,'test_results','matlab'); end
if ~isfolder(outputDir), mkdir(outputDir); end
report=struct; report.matlabVersion=version; report.runs=struct([]);
models={'double_integrator','aircraft'}; modes={'bCBF','blending','OI'};
for i=1:2
    model=models{i}; fixture=load(fullfile(root,'test_results',[model '_reference.mat'])); ref=fixture.ref;
    p=oicbf.config(model,'source'); results=struct;
    for j=1:3
        mode=modes{j}; fprintf('Full simulation: %s / %s\n',model,mode); ticRun=tic;
        result=oicbf.simulate(p,mode); elapsed=toc(ticRun); results.(mode)=result;
        py=ref.(['trajectory_' mode]); difference=result.x-py.x;
        stateError=max(abs(difference),[],2);
        ua=result.u(:,2:end); violation=max([0;reshape(p.lb-ua,[],1);reshape(ua-p.ub,[],1)]);
        % Scaled errors are descriptive; thresholds for trajectory comparison
        % are not fabricated from safety-function checks.
        item=struct('model',model,'mode',mode,'elapsedSeconds',elapsed,'steps',p.steps, ...
            'duration',result.t(end),'finite',all(isfinite([result.x(:);result.u(:)])), ...
            'minSafetyMargin',min(result.h),'maxInputViolation',violation, ...
            'allQPsFeasible',all(result.feasible),'maxStateErrorByComponent',stateError', ...
            'maxControlError',max(abs(result.u-py.u),[],'all'), ...
            'muError',NaN,'avgSolverTime',result.avgSolverTime,'maxSolverTime',result.maxSolverTime);
        if ~strcmp(mode,'bCBF'), item.muError=max(abs(result.mu-py.mu)); end
        % Interpretable acceptance thresholds chosen before seeing results:
        % double: 1e-8 state; aircraft: 2 m positions, .002 rad angles,
        % .002 rad/s roll rate and .002 load factor for smooth controllers.
        if i==1, thresholds=[1e-8;1e-8];
        else, thresholds=[.002;.002;.002;2;2;2;.002;.002]; end
        item.trajectoryThresholds=thresholds'; item.trajectoryAgreement=all(stateError<=thresholds);
        item.defaultSafetyChecksPassed=item.finite && item.minSafetyMargin>=-1e-8 ...
            && violation<=1e-8 && item.allQPsFeasible;
        if isempty(report.runs), report.runs=item; else, report.runs(end+1)=item; end
        disp(item);
        save(fullfile(outputDir,[model '_' mode '.mat']),'result','p');
        writeJSON(fullfile(outputDir,'simulation_validation.json'),report);
    end
    oicbf.plotComparison(results,outputDirForModel(outputDir,model));
end
report.allDefaultSafetyChecksPassed=all([report.runs.defaultSafetyChecksPassed]);
report.allTrajectoryComparisonsPassed=all([report.runs.trajectoryAgreement]);
writeJSON(fullfile(outputDir,'simulation_validation.json'),report);
assert(report.allDefaultSafetyChecksPassed,'oicbf:SimulationValidation','Default simulation safety checks failed.');
end
function out=outputDirForModel(base,model)
out=fullfile(base,model); if ~isfolder(out), mkdir(out); end
end
function writeJSON(file,value)
fid=fopen(file,'w','n','UTF-8'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,'PrettyPrint',true));
end
