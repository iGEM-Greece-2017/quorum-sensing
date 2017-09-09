function result= solveProblem(model,tlist,params)

model.SolverOptions.AbsoluteTolerance= params.solve.AbsTol;
model.SolverOptions.RelativeTolerance= params.solve.RelTol;
model.SolverOptions.ResidualTolerance= params.solve.RelTol;
model.SolverOptions.MaxIterations= 40;
model.SolverOptions.MinStep= params.solve.FeatureSize/10;
model.SolverOptions.ReportStatistics= 'on';

global bactRingDensity;
bactRingDensity= fewcell.util.bactRingDensity(params.g.bactCenters(:,1),params.g.bactSize, params.g.lateralSpacing);
totalBacteria= round(sum(bactRingDensity));
fprintf('Total bacteria: %d\n', totalBacteria);
% Defines the global variable <bactNodeCoeffs>, determining how the singlecell results 
%   affect the mesh nodes
fewcell.util.makeBactNodeCoeffs(model.Mesh,params.g.bactSize,bactRingDensity);

% pass params to <solveTimeDependent.m> without changing all the functions in between
global solveInternalParams;
solveInternalParams.y0= params.solve.y0;
solveInternalParams.AbsTol_y= params.solve.AbsTol_y;
tic;
result= solvepde(model,tlist);
toc;
global yyResults;
for b= 1:length(yyResults)
  yyResults{b}= [yyResults{b},totalBacteria*ones(length(tlist),1)];
end
