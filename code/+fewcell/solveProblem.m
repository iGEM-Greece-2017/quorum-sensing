function result= solveProblem(model,tlist,params)

model.SolverOptions.AbsoluteTolerance= params.AbsTol;
model.SolverOptions.RelativeTolerance= params.RelTol;
model.SolverOptions.ResidualTolerance= 1e-7;
model.SolverOptions.MaxIterations= 40;
model.SolverOptions.MinStep= params.FeatureSize/3;
model.SolverOptions.ReportStatistics= 'on';
tic;
result= solvepde(model,tlist);
toc;
