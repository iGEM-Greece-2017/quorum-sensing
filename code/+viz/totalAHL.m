function totalAHL= totalAHL(hist,winL)

totalAHL= arrayfun(@(x) sum(x.AHL(:)),hist);
totalAHL(:,2)= util.zerophMAfilt(winL,totalAHL);
plot(totalAHL);
grid minor; title('Total extracellular [AHL] in grid');
