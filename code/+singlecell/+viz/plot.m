function plot(t,y, figID,logplot)

t_h= t/60;
if logplot      % cutoff log plots at a minimum value
  logCutoff= 10.^(log10(max(y(:,[1,4,5,8,9]),[],1)) - [6,4,4,6,6]);
  y(:,[1,4,5,8,9])= max(y(:,[1,4,5,8,9]), logCutoff);
end

figure(figID); clf;
subplot(221);
if ~logplot, yyaxis left; plot(t_h, y(:,1)); ylabel('c [nM]'); axis tight;
else, yyaxis left; semilogy(t_h, y(:,1)); ylabel('c [nM]'); axis tight;
end
yyaxis right; semilogy(t_h, y(:,11)); ylabel('cells'); axis tight;
xlabel('t [hour]'); xlim([t_h(1) t_h(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('DNA'); legend('DNA', 'cell number','location','southeast');

subplot(222); yyaxis left; plot(t_h, y(:,6)); ylabel('c [nM]'); axis tight;
yyaxis right; plot(t_h, y(:,10)); ylabel('c [nM]'); axis tight;
xlabel('t [hour]'); xlim([t_h(1) t_h(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('AHL'); legend('AHL_{cyt}', 'AHL_{agar}', 'location','southeast');

subplot(223);
if ~logplot
  yyaxis left; plot(t_h,y(:,5)); ylabel('c [nM]'); axis tight;
  yyaxis right; plot(t_h,y(:,4)); ylabel('c [nM]'); axis tight;
else
  yyaxis left; semilogy(t_h,y(:,5)); ylabel('c [nM]'); axis tight;
  yyaxis right; semilogy(t_h,y(:,4)); ylabel('c [nM]'); axis tight;
end
xlabel('t [hour]'); xlim([t_h(1) t_h(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('LuxR, LuxI'); legend('LuxR','LuxI', 'location','southeast');

subplot(224);
if ~logplot
  yyaxis left; plot(t_h,y(:,8)); ylabel('c [nM]'); axis tight;
  yyaxis right; plot(t_h,y(:,9)); ylabel('c [nM]'); axis tight;
else
  yyaxis left; semilogy(t_h,y(:,8)); ylabel('c [nM]'); axis tight;
  yyaxis right; semilogy(t_h,y(:,9)); ylabel('c [nM]'); axis tight;
end
xlabel('t [hour]'); xlim([t_h(1) t_h(end)]);
g= gca; g.XMinorGrid= 'on'; g.YMinorGrid= 'on';
title('LuxRAHL2, DNALuxRAHL2'); legend('LuxRAHL2','DNALuxRAHL2', 'location','southeast');
