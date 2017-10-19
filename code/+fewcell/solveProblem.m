function result= solveProblem(model,tlist,params)

model.SolverOptions.AbsoluteTolerance= params.solve.AbsTol;
model.SolverOptions.RelativeTolerance= params.solve.RelTol;
model.SolverOptions.ResidualTolerance= params.solve.RelTol;
model.SolverOptions.MaxIterations= 40;
model.SolverOptions.MinStep= params.solve.FeatureSize/20;
model.SolverOptions.ReportStatistics= 'off';

% pass params to <solveTimeDependent.m> without changing all the functions in between
global solveInternalParams;
nBact= size(params.g.bactCenters,1);
if length(params.solve.y0)==8
  solveInternalParams.y0= repmat(params.solve.y0,nBact,1);
else
  solveInternalParams.y0= params.solve.y0;
end
solveInternalParams.AbsTol_y= params.solve.AbsTol_y;
result= solvepde(model,tlist);
global yyResults;
global bactNodes;
for b= 1:length(yyResults)
  bNode= find(bactNodes(:,b),1);
  cAHL= result.NodalSolution([bNode,1],:)';   % Assume: 1: agar surface outside, 6: bact1
  yyResults{b}(:,[6,10])= cAHL;
  yyResults{b}= [yyResults{b},ones(length(tlist),1)*nBact];
end
