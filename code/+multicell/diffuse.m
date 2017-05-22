function AHL= diffuse(AHL, diffusionConstant)
% Diffuse AHL in a hexagonal neighborhood according to a 2-point diffusion equation. The base diffusion
% equation is 1D and unidirectional, representing the diffusion between 2 neighbors in 1 direction.
% However each cell has 6 neighbors arranged in a circular topology (simulating hex cells),
% so rotational symmetry is then applied and that same amount goes to each neighbor. In case more
% AHL is to be sent out in total than is available, the amount is truncated.
% 
% The calculation occurs in 2 steps to simulate the full bidirectional diffusion.
% 1. Each cell, given only its own AHL, calculates the total amount to give to all neighbors,
%    which is subtracted from its own amount. Thus, 1 side of the diffusion equation is implemented.
% 2. AHL is traded: each cell gives 1/6 of its amount-to-give to each neighbor, thus implementing
%    the reciprocal side of the diffusion equation
% 
% Note: GPU-friendly
% Warning: this algorithm will become unstable for large timesteps

%% Step 1: calculate amounts-to-give
% toGive [Nx]x[Ny]: total amount to give to all neighbors
toGive= % TODO

%% Step 2: filter image to add neighbors' toGive values and to remove own value
% Actually, because a hexagonal grid is simulated, each row must be alternately treated with
% a different filter shape
diffusionFilter.left= [1 1 0; 1 -6 1; 1 1 0]/6;
diffusionFilter.right= [0 1 1; 1 -6 1; 0 1 1]/6;
AHL_left= imfilter(AHL,diffusionFilter.left);
AHL_right= imfilter(AHL,diffusionFilter.right);
% Fill alternate rows with the previous results
AHL(1:2:end,:)= AHL_left(1:2:end,:);
AHL(2:2:end,:)= AHL_right(2:2:end,:);
end
