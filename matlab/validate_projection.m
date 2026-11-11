function report = validate_projection()
%VALIDATE_PROJECTION Compare custom projection against independent quadprog.
root=fileparts(fileparts(mfilename('fullpath')));
out=fullfile(root,'test_results','matlab'); if ~isfolder(out), mkdir(out); end
data=load(fullfile(root,'test_results','projection_reference.mat'));
report=struct('cases',numel(data.cases),'feasibilityMismatches',0,'maxControlError',0, ...
    'maxConstraintViolation',0,'passed',false);
for k=1:numel(data.cases)
    c=data.cases{k}; n=double(c.n);
    [u,ok]=oicbf.projectPolytope(c.up(:),c.C,c.d(:),-ones(n,1),ones(n,1));
    report.feasibilityMismatches=report.feasibilityMismatches+(ok~=c.feasible);
    if ok && c.feasible
        report.maxControlError=max(report.maxControlError,max(abs(u-c.u(:))));
        report.maxConstraintViolation=max(report.maxConstraintViolation,max(c.d(:)-c.C*u));
    end
end
report.passed=report.feasibilityMismatches==0 && report.maxControlError<1e-8 && report.maxConstraintViolation<1e-8;
fid=fopen(fullfile(out,'projection_validation.json'),'w');
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,'%s\n',jsonencode(report,'PrettyPrint',true));
disp(report); assert(report.passed,'oicbf:ProjectionValidation','Projection differs from quadprog oracle.');
end
