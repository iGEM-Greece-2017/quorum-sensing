function [bacteriumState, AHL]= singleBactModel(hasBacterium,bacteriumState,AHL, tstep, params)
% Applies a single bacterium model to each bacterium on a cell grid individually
% Input:
% - hasBacterium [Nx]x[Ny] {sparse<->full,logical}
% - bacteriumState [Nx]x[Ny] {sparse<->full}
% - AHL [Nx]x[Ny]: AHL concentration. In every cell with a bacterium, it corresponds to 
%     single cell model's outside AHL concentration

t_s= tstep*params.dt; t_f= (tstep+1)*params.dt;
% Execute single bacterium model foreach grid cell
[bacteriumState(hasBacterium), AHL(hasBacterium)]= ...
  arrayfun(@(bs,ahl) runSingleBact(bs,ahl,t_s,t_f), ...
           bacteriumState(hasBacterium), AHL(hasBacterium));
end

function [bs, ahl]= runSingleBact(bs,ahl, t_s,t_f)
t = linspace(t_s/60,t_f/60, ceil((t_f-t_s)*100)+3);  % Times from seconds to minutes
y= ode4(@singlecell.model,t,[bs{:};ahl]);
bs= {y(end,1:9)'}; ahl= y(end,10);
end
