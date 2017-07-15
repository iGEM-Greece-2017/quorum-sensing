function plot(t,y0,y,N, figID)
t= t/3600;
%titles= sprintf(', DNA=%.2gM, y_6(0)=%.2gM', y0(1),y0(6));
figure(figID);
subplot(221); plot(t,y(:,1)); grid minor; axis tight; 
title('DNA'); ylabel('c [M]'); xlabel('t [hour]');

subplot(222); AHL= y(:,6);
plot(t,AHL); grid minor; axis tight; 
title('AHL');
%title(['AHL, ',num2str(N),' cell',titles]);
ylabel('c [M]'); xlabel('t [hour]');

subplot(223); yyaxis left; plot(t,y(:,5)); axis tight; ylabel('c [M]');
yyaxis right; plot(t,y(:,4)); axis tight; ylabel('c [M]'); xlabel('t [hour]');
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('LuxR, LuxI');
legend('LuxR','LuxI', 'location','southeast');

subplot(224); yyaxis left; plot(t,y(:,8)); axis tight; ylabel('c [M]');
yyaxis right; plot(t,y(:,9)); axis tight; ylabel('c [M]'); xlabel('t [hour]');
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['LuxRAHL2, DNALuxRAHL2']);
legend('LuxRAHL2','DNALuxRAHL2', 'location','southeast');

suptitle(sprintf('N=%.1e cells, init DNA=%.2gM', N,y0(1)));
