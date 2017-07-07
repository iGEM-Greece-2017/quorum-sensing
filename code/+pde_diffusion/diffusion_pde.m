% Setup the diffusion pde
import pde_diffusion.*;
%% Geometry definition
[bact(:,1),name{1}]= bactGeometry([-3.1,-4.3],0.2,1);
[bact(:,2),name{2}]= bactGeometry([-4.7,-3.9],0.2,2);
[bact(:,3),name{3}]= bactGeometry([-3.2, 2.1],0.2,3);
[bact(:,4),name{4}]= bactGeometry([ 2.5,-0.6],0.2,4);

domainLim= 30;
domain= [3;4;-domainLim;-domainLim;domainLim;domainLim;-domainLim;domainLim;domainLim;-domainLim];

geomNames= char('domain',name{:})';
setFormula= [geomNames(:,1)','-',geomNames(:,2)','-',geomNames(:,3)','-',geomNames(:,4)','-',geomNames(:,5)'];
%setFormula= ['domain-',name{1},'-',name{2},'-',name{3},'-',name{4}];
dl= decsg([domain,bact], setFormula, geomNames);

%% PDE
numberOfPDE= 1;
model= createpde(numberOfPDE);
geometryFromEdges(model,dl);

figure(1);
pdegplot(model,'EdgeLabels','on','FaceLabels','on');
drawnow;
axis equal

%% Boundary conditions
applyBoundaryCondition(model,'neumann','Edge',[1,2,11,12],'g',0,'q',1);
applyBoundaryCondition(model,'neumann','Edge',[3:10,13:20],'g',1,'q',1e-3);

%% PDE coeffs
specifyCoefficients(model, 'm',0, 'd',1, 'c',1, 'a',0, 'f',0);  % d: d*u', c: -div(c*grad(u)), f: source

%% Initial conditions
setInitialConditions(model,0);

%% Mesh
msh= generateMesh(model,'Hgrad',1.2,'Hmax',1.8);
figure(2);
pdemesh(model);
drawnow;
axis equal

%% Time discretization
nframes = 40;
tlist = logspace(-1,2,nframes);

%% Solve
model.SolverOptions.ReportStatistics ='on';
result= solvepde(model,tlist);
u1 = result.NodalSolution;

%% Plot solution
figure(3);
umax = max(max(u1));
umin = min(min(u1));
for j = 1:nframes
  pdeplot(model,'XYData',u1(:,j),'ZData',u1(:,j));
  caxis([umin umax]);
  xlim([-10 10]); ylim([-10 10]);
  drawnow;
  pause(0.05);
end
