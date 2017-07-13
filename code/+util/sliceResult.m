function r= sliceResult(model,result,t)
% Hack Matlab and take a slice out of the TimeDependentResults object as a StationaryResults

r= StationaryResults(model, result.NodalSolution(:,1), {result,t}, 'hackit, please!');
