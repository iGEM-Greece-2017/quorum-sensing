function [Q,G,H,R]=assembleBoundary(self,varargin)
%assembleBoundary Assembles boundary condition contributions in a PDE problem.
%

%       Copyright 2015-2016 The MathWorks, Inc.

narginchk(1,3);
narg = nargin;

gotu=false;
gottime=false;

msh = self.Mesh;
[p,e,~] = msh.meshToPet();
bl = pde.PDEModel(self.PDESystemSize);
if isempty(self.BoundaryConditions)
    bc = [];  
else
    bc = self.BoundaryConditions.packBoundaryConditions();
end

if narg>1
    u = varargin{1};
    gotu=true;
end

if narg == 3
    time=varargin{2};
    gottime=true;
end

if ~gotu
    u = [];
end

if ~gottime
    time = [];
end

npr = size(p,1);
if(npr == 3)
    bcImpl = pde.internal.pde3DBCImpl(bl,bc,true); % true for R2016b BC implementation
    [Q,G,H,R] = bcImpl.getBCMatrices(p,e,u,time);
    return;
elseif(npr == 2)
    bcImpl = pde.internal.pde2DBCSfnImpl(bl,bc,true); % true for R2016b BC implementation
    [Q,G,H,R] = bcImpl.getBCMatrices(p,e,u,time,msh.GeometricOrder);
    return
else
    error(message('pde:pdeModel:invalidP', npr));
end

end


