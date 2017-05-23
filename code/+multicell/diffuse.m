function AHL= diffuse(AHL,dt,r, const)
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
toGive= AHL*const.t_s ./ (2*dt) .* exp(-(2*r)^2*const.t_s*pi / (2*dt));

%% Step 2: filter image to add neighbors' toGive values and to remove own value
% Actually, because a hexagonal grid is simulated, each row must be alternately treated with
% a different filter shape
diffusionFilter.left= [1 1 0; 1 -6 1; 1 1 0]/6;
diffusionFilter.right= [0 1 1; 1 -6 1; 0 1 1]/6;
dAHL_left= imfilter(toGive,diffusionFilter.left);      % Only correct for half the rows
dAHL_right= imfilter(toGive,diffusionFilter.right);    % Only correct for the other half
% Fill alternate rows with the previous results
AHL(1:2:end,:)= AHL(1:2:end,:) + dAHL_left(1:2:end,:);
AHL(2:2:end,:)= AHL(2:2:end,:) + dAHL_right(2:2:end,:);
end
