function [q,tstep]= quantizeGC(x,q0,qsteps)
% Calculates the time points where the integration should be interrupted to accomodate
% bacterial population growth

qLevels= cumsum([q0,qsteps]);
q= x;
qi= 1;
for t= 1:length(x)
  if qi < length(qLevels)
    qi= qi + (abs(x(t)-qLevels(qi)) > abs(x(t)-qLevels(qi+1))); 
  end
  q(t)= qLevels(qi);
end
tstep= find(diff(q));
end
