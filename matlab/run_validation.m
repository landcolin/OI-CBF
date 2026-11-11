function report = run_validation()
%RUN_VALIDATION Run formula checks and full simulations with original fixtures.
report.math=validate_math();
report.projection=validate_projection();
report.paper=validate_paper_profile();
report.simulations=validate_simulations();
if ~isempty(which('quadprog')), report.toolbox=validate_toolbox(); end
if ~isempty(which('sym')), report.symbolic=validate_symbolic(); end
if ~isempty(which('parpool')), report.parallel=validate_parallel(); end
end
