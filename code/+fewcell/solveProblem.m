function result= solveProblem(model,tlist,params)

model.SolverOptions.AbsoluteTolerance= params.AbsTol;
model.SolverOptions.RelativeTolerance= params.RelTol;
%model.SolverOptions.ResidualTolerance= params.RelTol;
model.SolverOptions.MaxIterations= 30;
model.SolverOptions.MinStep= params.FeatureSize/10;
model.SolverOptions.ReportStatistics= 'on';
% pass params to <solveTimeDependent.m> without changing all the functions in between
global solveInternalParams;
solveInternalParams.y0= params.y0;
solveInternalParams.AbsTol_y= params.AbsTol_y;
tic;
result= solvepde(model,tlist);
toc;
