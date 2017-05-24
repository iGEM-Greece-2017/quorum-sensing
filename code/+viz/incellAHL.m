function incellAHL= incellAHL(hist,cellnum)

cells= hist(1).hasBact; 
sel= find(cells); sel= sel(floor(rand(cellnum,1)*(length(sel)-0.01))+1);
unsel= setdiff(find(cells),sel);
cells(unsel)= 0;
%% Gather intracellular AHL
% incellAHL: [t] cellarray with [sum(cells)] matrices
incellAHL= arrayfun(@(x) cellfun(@(x) x(6), x.bactState(cells)), hist, 'UniformOutput',false);
incellAHL= cell2mat(incellAHL')'; % [t]x[sum(cells)]

%% Plot
plot(incellAHL);
grid minor; title('Intracellular AHL for some bacteria');
end
