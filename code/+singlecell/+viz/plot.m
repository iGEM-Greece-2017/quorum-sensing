function plot(t,y,N, figID)
t_min= t/60;
%titles= sprintf(', DNA=%.2gM, y_6(0)=%.2gM', y0(1),y0(6));
figure(figID);
subplot(221); plot(t_min,y(:,1)./N); grid minor; axis tight;
title('DNA (per cell)'); ylabel('c [nM]'); xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);

subplot(222); yyaxis left; plot(t_min, y(:,6)); ylabel('c [nM]'); axis tight;
yyaxis right; plot(t_min, y(:,10)); ylabel('c [nM]'); axis tight;
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('AHL'); legend('AHL_{cyt}', 'AHL_{agar}', 'location','southeast');


subplot(223); yyaxis left; plot(t_min,y(:,5)./N); ylabel('c [nM]'); axis tight;
yyaxis right; plot(t_min,y(:,4)./N); ylabel('c [nM]'); axis tight;
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('LuxR, LuxI (per cell)'); legend('LuxR','LuxI', 'location','southeast');

subplot(224); yyaxis left; plot(t_min,y(:,8)./N); ylabel('c [nM]'); axis tight;
yyaxis right; plot(t_min,y(:,9)./N); ylabel('c [nM]'); axis tight;
xlabel('t [hour]'); xlim([t_min(1) t_min(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['LuxRAHL2, DNALuxRAHL2 (per cell)']); legend('LuxRAHL2','DNALuxRAHL2', 'location','southeast');

suptitle(sprintf('N=%d cells', N));
