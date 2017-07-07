%% Runs the single cell model
t= linspace(0,1000,10000);
y0= [1.6E-7;0;0;0;0;2E-5;0;0;0;0];
opt= odeset('Jacobian', @singlecell.modelJacobian, 'RelTol',1e-5,'AbsTol',1e-10,'InitialStep',1e-6);
[t,y]= ode23tb(@singlecell.model,t,y0,opt);

figure(1);
subplot(221); plot(t,y(:,1)); grid minor; title('y_1 for DNA=1.6E-7, y_6(0)=2E-5'); subplot(222); plot(t,y(:,6)); grid minor; title('y_6 for DNA=1.6E-7, y_6(0)=2E-5');
subplot(223); plot(t,y(:,2)); grid minor; title('y_2 for DNA=1.6E-7, y_6(0)=2E-5'); subplot(224); plot(t,y(:,8:9)); grid minor; title('y_8,y_9 for DNA=1.6E-7, y_6(0)=2E-5'); legend('y_8','y_9');
