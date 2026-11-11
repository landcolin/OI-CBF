function results = run_comparison(model,profile,outputDir,qpSolver)
%RUN_COMPARISON Run all three methods and generate native MATLAB figures.
% addpath('matlab'); run_comparison('double_integrator');
% run_comparison('aircraft','source');
if nargin<1, model='double_integrator'; end
if nargin<2, profile='source'; end
if nargin<3 || isempty(outputDir), outputDir=fullfile(fileparts(mfilename('fullpath')),'results',model,profile); end
if nargin<4, qpSolver='projection'; end
if ~isfolder(outputDir), mkdir(outputDir); end
p=oicbf.config(model,profile); p.qpSolver=validatestring(qpSolver,{'projection','quadprog'});
modes={'bCBF','blending','OI'}; results=struct;
for i=1:numel(modes)
    mode=modes{i}; fprintf('Running %s / %s / %s\n',model,profile,mode);
    results.(mode)=oicbf.simulate(p,mode);
end
save(fullfile(outputDir,'comparison.mat'),'results','p');
oicbf.plotComparison(results,outputDir);
end
