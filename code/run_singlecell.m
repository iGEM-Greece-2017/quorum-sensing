% For unknown reasons, this doesn't work well
%[t,y]= ode15s(@singlecell.model,[0 1000],[1.6E-9;0;0;0;0;2E-5;0;0;0;0]);

% Finally runs well:
t = linspace(0,800,80000);
y= ode4(@singlecell.model,t,[1.6E-7;0;0;0;0;2E-5;0;0;0;0]);

subplot(221); plot(t,y(:,1)); grid minor; title('y_1 for DNA=1.6E-7, y_6(0)=2E-5'); subplot(222); plot(t,y(:,6)); grid minor; title('y_6 for DNA=1.6E-7, y_6(0)=2E-5');
subplot(223); plot(t,y(:,2)); grid minor; title('y_2 for DNA=1.6E-7, y_6(0)=2E-5'); subplot(224); plot(t,y(:,8:9)); grid minor; title('y_8,y_9 for DNA=1.6E-7, y_6(0)=2E-5'); legend('y_8','y_9');
