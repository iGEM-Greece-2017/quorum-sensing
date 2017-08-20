function plot(t,y,N0, figID,perCell)

cellNorm= ones(length(t),1);
cellnormTitle= ' (total)';
if perCell
  cellNorm= y(:,11);
  cellnormTitle= ' (per cell)';
end
  
t_min= t/60;
%titles= sprintf(', DNA=%.2gM, y_6(0)=%.2gM', y0(1),y0(6));
figure(figID);
subplot(221);
if perCell, yyaxis left; plot(t_min, y(:,1)./cellNorm); ylabel('c [nM]'); axis tight;
else, yyaxis left; semilogy(t_min, y(:,1)./cellNorm); ylabel('c [nM]'); axis tight;
end
yyaxis right; semilogy(t_min, y(:,11)); ylabel('cells'); axis tight;
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['DNA',cellnormTitle]); legend('DNA', 'cell number','location','southeast');

subplot(222); yyaxis left; plot(t_min, y(:,6)); ylabel('c [nM]'); axis tight;
yyaxis right; plot(t_min, y(:,10)); ylabel('c [nM]'); axis tight;
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('AHL'); legend('AHL_{cyt}', 'AHL_{agar}', 'location','southeast');

subplot(223);
if perCell
  yyaxis left; plot(t_min,y(:,5)./cellNorm); ylabel('c [nM]'); axis tight;
  yyaxis right; plot(t_min,y(:,4)./cellNorm); ylabel('c [nM]'); axis tight;
else
  yyaxis left; semilogy(t_min,y(:,5)./cellNorm); ylabel('c [nM]'); axis tight;
  yyaxis right; semilogy(t_min,y(:,4)./cellNorm); ylabel('c [nM]'); axis tight;
end
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['LuxR, LuxI',cellnormTitle]); legend('LuxR','LuxI', 'location','southeast');

subplot(224);
if perCell
  yyaxis left; plot(t_min,y(:,8)./cellNorm); ylabel('c [nM]'); axis tight;
  yyaxis right; plot(t_min,y(:,9)./cellNorm); ylabel('c [nM]'); axis tight;
else
  yyaxis left; semilogy(t_min,y(:,8)./cellNorm); ylabel('c [nM]'); axis tight;
  yyaxis right; semilogy(t_min,y(:,9)./cellNorm); ylabel('c [nM]'); axis tight;
end
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['LuxRAHL2, DNALuxRAHL2',cellnormTitle]); legend('LuxRAHL2','DNALuxRAHL2', 'location','southeast');

%suptitle(sprintf('N0=%d cells', N0));
