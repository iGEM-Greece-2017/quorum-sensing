function [tstep,tLast]= makeGrowthTimesteps(growth)
% Calculates the time points where the integration should be interrupted to accomodate
% bacterial population growth
  tstep= [];
  i= 1;
  tPrev= 1;
  tLast= length(growth.bactN);
  grownThisStep= false;
  for t= 2:length(growth.bactN)
    if growth.bactN(t)-growth.bactN(tPrev) > growth.stepSize(i)/3   % Timestamp the middle of the growth step, not the end
      if ~grownThisStep, tstep= [tstep,t]; end
      grownThisStep= true;
      if growth.bactN(t)-growth.bactN(tPrev) > growth.stepSize(i)
        i= i+1;
        tPrev= t;
        grownThisStep= false;
        if i > length(growth.stepSize) && t<length(growth.bactN)
          warning('[makeGrowthTimesteps]: Geometry too small to accomodate all the bacteria!');
          tLast= t;
          break;
        end
      end
    end
  end
  if isempty(tstep)
    error('[makeGrowthTimesteps]: No growth timesteps! Number of rings is probably too small or growth too slow');
  end
  if tstep(end) < length(growth.bactN)
    tstep= [tstep,length(growth.bactN)];
  end
end
