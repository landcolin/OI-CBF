function report = validate_parallel()
%VALIDATE_PARALLEL Verify serial/parallel reproducibility including failures.
serial=run_batch_tests(3,'serial'); parallel=run_batch_tests(3,'parallel',2);
a=serial.results; b=parallel.results; mismatch=0; maxDifference=0;
for k=1:numel(a)
    same=strcmp(a(k).status,b(k).status) && strcmp(a(k).errorId,b(k).errorId);
    mismatch=mismatch+~same;
    if strcmp(a(k).status,'completed') && strcmp(b(k).status,'completed')
        maxDifference=max(maxDifference,max(abs(a(k).finalState-b(k).finalState)));
        maxDifference=max(maxDifference,abs(a(k).minSafetyMargin-b(k).minSafetyMargin));
        maxDifference=max(maxDifference,abs(a(k).maxInputViolation-b(k).maxInputViolation));
        mismatch=mismatch+(a(k).finite~=b(k).finite);
    end
end
report=struct('cases',numel(a),'outcomeMismatches',mismatch,'maxNumericalDifference',maxDifference, ...
    'serialSeconds',serial.computeSeconds,'parallelSeconds',parallel.computeSeconds, ...
    'poolStartupSeconds',parallel.poolStartupSeconds,'workers',parallel.workers, ...
    'completed',parallel.completed,'infeasible',parallel.infeasible,'otherErrors',parallel.otherErrors, ...
    'safetyViolations',parallel.safetyViolations,'inputViolations',parallel.inputViolations);
completed=strcmp({b.status},'completed');
report.passed=mismatch==0 && maxDifference<1e-9 && serial.otherErrors==0 && parallel.otherErrors==0 ...
    && all([b(completed).finite]);
root=fileparts(fileparts(mfilename('fullpath')));
fid=fopen(fullfile(root,'test_results','matlab','parallel_validation.json'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true)); disp(report);
assert(report.passed,'oicbf:ParallelValidation','Serial and parallel batch results differ.');
end
