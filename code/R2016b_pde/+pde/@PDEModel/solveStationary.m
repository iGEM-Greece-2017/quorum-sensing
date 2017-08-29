function u = solveStationary(self,coefstruct,varargin)
% solveStationary - Internal function for solving a time-independent PDE.

if isempty(varargin)
    % Constant coefficients case, compute FE matrices.
    [p,~,t] = self.Mesh.meshToPet();
    thePde = pde.PDEModel(self.PDESystemSize);
    thePde.BoundaryConditions = self.BoundaryConditions;
    ndims = size(p,1);
    % Caculate Global FE matrices
    if(ndims==2)
        [K, F] = formGlobalKF2D(thePde, p, t, coefstruct,[],[]);
        A = formGlobalM2D(thePde, p, t, coefstruct,[],[],'a');
    elseif(ndims==3)
        [K, F] = formGlobalKF3D(thePde, p, t, coefstruct,[],[]);
        A = formGlobalM3D(thePde, p, t, coefstruct,[],[],'a');
    end
    [Q,G,H,R]=self.assembleBoundary();
else
    %Coefficients might depended on spatial coordinates, reuse the system
    %matrices available in DiscretizedPDEModel object passed in as varargin.
    femodel = varargin{1};
    K = femodel.K;
    A = femodel.A;
    F = femodel.F;
    Q = femodel.Q;
    G = femodel.G;
    H = femodel.H;
    R = femodel.R;
end

[null,orth]=pdenullorth(H);
if size(orth,2)==0
    ud=zeros(size(K,2),1);
else
    ud=full(orth*((H*orth)\R));
end

KK=K+A+Q;
FF=null'*((F+G)-KK*ud);
KK=null'*KK*null;
KK = self.checkForSymmetry(KK);

if size(null,2)==0
    u=ud;
else
    u=null*(KK\FF)+ud;
end

end
