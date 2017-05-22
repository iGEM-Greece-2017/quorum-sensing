function [bacteriumState, AHL]= singleBactModel(hasBacterium,bacteriumState,AHL, params)
% Applies a single bacterium model to each bacterium on a cell grid individually
% Input:
% - hasBacterium [Nx]x[Ny] {sparse<->full,logical}
% - bacteriumState [Nx]x[Ny] {sparse<->full}
% - AHL [Nx]x[Ny]: AHL concentration. In every cell with a bacterium, it corresponds to 
%     single cell model's outside AHL concentration
% - params: single bacterium model parameters

[bacteriumState(hasBacterium), AHL(hasBacterium)]= arrayfun(@(state,AHL) singleBact(state,AHL,params),...
    bacteriumState(hasBacterium), AHL(hasBacterium));
end
