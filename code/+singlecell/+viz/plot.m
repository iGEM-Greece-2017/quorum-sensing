function plot(t,y0,y, figID)
titles= sprintf(' for DNA=%.2gM, y_6(0)=%.2gM', y0(1),y0(6));
figure(figID);
subplot(221); plot(t,y(:,1)); grid minor; axis tight; 
title(['y_1',titles]); ylabel('c [M]'); xlabel('t [s]');
subplot(222); plot(t,y(:,6)); grid minor; axis tight; 
title(['y_6',titles]); ylabel('c [M]'); xlabel('t [s]');
subplot(223); plot(t,y(:,2)); grid minor; axis tight; 
title(['y_2',titles]); ylabel('c [M]'); xlabel('t [s]');
subplot(224); yyaxis left; plot(t,y(:,8)); axis tight; ylabel('c [M]'); xlabel('t [s]');
yyaxis right; plot(t,y(:,9)); axis tight; ylabel('c [M]'); xlabel('t [s]');
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title(['y_8,y_9',titles]); legend('y_8','y_9');
