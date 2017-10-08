function [position,isterminal,direction]= events(t,~)
  global growth;
  position= t-growth.tstep(growth.i);
  isterminal= 1;
  direction= 0;
end
