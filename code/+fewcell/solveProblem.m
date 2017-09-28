function result= solveProblem(model,tlist,params)

model.SolverOptions.AbsoluteTolerance= params.solve.AbsTol;
model.SolverOptions.RelativeTolerance= params.solve.RelTol;
model.SolverOptions.ResidualTolerance= params.solve.RelTol;
model.SolverOptions.MaxIterations= 40;
model.SolverOptions.MinStep= params.solve.FeatureSize/20;
model.SolverOptions.ReportStatistics= 'on';

%    ~~~ NOT ANYMORE ~~~
% Defines the global variable <bactNodeCoeffs>, determining how the singlecell results 
%   affect the mesh nodes
nBact= size(params.g.bactCenters,1);
fewcell.util.makeBactNodeCoeffs(model.Mesh,nBact);

% pass params to <solveTimeDependent.m> without changing all the functions in between
global solveInternalParams;
solveInternalParams.y0= params.solve.y0;
solveInternalParams.AbsTol_y= params.solve.AbsTol_y;
tic;
result= solvepde(model,tlist);
toc;
cAHL= result.NodalSolution([6,1],:)';   % Assume: 1: agar surface outside, 6: bact1
global yyResults;
for b= 1:length(yyResults)
  yyResults{b}(:,[6,10])= cAHL;
  yyResults{b}= [yyResults{b},ones(length(tlist),1)];
end
