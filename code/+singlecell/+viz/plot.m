function plot(t,y,N, figID)
t= t/3600;
%titles= sprintf(', DNA=%.2gM, y_6(0)=%.2gM', y0(1),y0(6));
figure(figID);
subplot(221); plot(t,y(:,1)./N); grid minor; axis tight;
title('DNA'); ylabel('c [M]'); xlabel('t [hour]'); xlim([t(1) t(end)]);

subplot(222); AHL= y(:,6);
plot(t,AHL); grid minor; axis tight; 
title('total AHL');
%title(['AHL, ',num2str(N),' cell',titles]);
ylabel('c [M]'); xlabel('t [hour]'); xlim([t(1) t(end)]);

subplot(223); yyaxis left; plot(t,y(:,5)./N); ylabel('c [M]'); axis tight;
yyaxis right; plot(t,y(:,4)); ylabel('c [M]'); axis tight;
xlabel('t [hour]'); xlim([t(1) t(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('LuxR, LuxI'); legend('LuxR','LuxI', 'location','southeast');

subplot(224); yyaxis left; plot(t,y(:,8)./N); ylabel('c [M]'); axis tight;
yyaxis right; plot(t,y(:,9)); ylabel('c [M]'); axis tight;
xlabel('t [hour]'); xlim([t(1) t(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['LuxRAHL2, DNALuxRAHL2']); legend('LuxRAHL2','DNALuxRAHL2', 'location','southeast');

suptitle(sprintf('N=%d cells', N));
