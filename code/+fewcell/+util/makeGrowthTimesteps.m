function [tstep,tLast]= makeGrowthTimesteps(growth)
% Calculates the time points where the integration should be interrupted to accomodate
% bacterial population growth
  tstep= [];
  i= 1;
  tPrev= 1;
  tLast= length(growth.bactN);
  for t= 2:length(growth.bactN)
    if growth.bactN(t)-growth.bactN(tPrev) > growth.stepSize(i)
      tstep= [tstep,t];
      tPrev= t;
      i= i+1;
      if i > length(growth.stepSize)
        warning('Geometry too small to accomodate all the bacteria!');
        tLast= t;
        break;
      end
    end
  end
end
