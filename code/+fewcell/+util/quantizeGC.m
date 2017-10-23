function [q,tstep,adapt_dr]= quantizeGC(x,q0,qsteps,min_dt0)
% Calculates the time points where the integration should be interrupted to accomodate
% bacterial population growth
% <adapt_dr> is the number of initial time steps concatenated into 1

  % Calculate initial quantized growth curve
  q= calcQGC(x,q0,qsteps);
  tstep= find(diff(q));
  tstep= [tstep;ones(length(qsteps)-length(tstep),1)*tstep(end)];
  adapt_dr= ones(size(tstep));
  dTstep= diff([0;tstep]);
  min_dt= min(max(min_dt0,2),tstep(1));
  
  % Concatenate multiple growth steps, if they are too small
  attempts= 0;
  while min(dTstep) < min_dt && attempts<100  % if >100 attempts, give up
    t= 1; t1= 1; concatFlag= false;
    while t<length(tstep)
      if ~concatFlag && dTstep(t) < min_dt
        t1= t-1;
        concatFlag= true;
      end
      if concatFlag && tstep(t) - tstep(t1) > min_dt
        qsteps= [qsteps(1:t1-1),sum(qsteps(t1:t-1)),qsteps(t:end)];
        adapt_dr= [adapt_dr(1:t1-1);sum(adapt_dr(t1:t-1));adapt_dr(t:end)];
        tstep= [tstep(1:t1);tstep(t:end)];
        dTstep= diff([0;tstep]);
        t= t1;
        concatFlag= false;
      end
      t= t+1;
    end
  
    % Calculate new quantized growth curve
    q= calcQGC(x,q0,qsteps);
    %tstep= find(diff(q));
    %dTstep= diff([0;tstep]);
    min_dt= min(max(min_dt0,2),tstep(1));
    attempts= attempts+1;
  end
  if any(dTstep<2 & ~diff([0;(dTstep<2)]))   % if the step time drops to 1 for consecutive timesteps...
    warning('[quantizeGC]: dr too small to follow the growth curve! Consider increasing it');
  end
  
  %fprintf('Quantization mre: %.2f%%\n', mean((x-q)./x)*100);
end

function q= calcQGC(x,q0,qsteps)
  qLevels= cumsum([q0,qsteps]);
  q= x;
  qi= 1;
  for t= 1:length(x)
    if qi < length(qLevels)
      qi= qi + (abs(x(t)-qLevels(qi)) > abs(x(t)-qLevels(qi+1))); 
    end
    q(t)= qLevels(qi);
  end
end
