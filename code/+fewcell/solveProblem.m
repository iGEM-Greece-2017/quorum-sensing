function result= solveProblem(model,tlist,params)

model.SolverOptions.AbsoluteTolerance= params.AbsTol;
model.SolverOptions.RelativeTolerance= params.RelTol;
%model.SolverOptions.ResidualTolerance= params.RelTol;
model.SolverOptions.MaxIterations= 30;
model.SolverOptions.MinStep= params.FeatureSize/10;
model.SolverOptions.ReportStatistics= 'on';
tic;
result= solvepde(model,tlist);
toc;
