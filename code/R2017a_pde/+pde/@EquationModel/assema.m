function [ks,km,fm] = assema(thePde, varargin) % thePde, p, t, c, a, f, varargin
%ASSEMA Calculate global matrices and load vector for 2D and 3D geometry
% Select between 2D and 3D implementations based on number of rows in p
%
% This undocumented function may be changed or removed in a future release.

%       Copyright 2014-2016 The MathWorks, Inc.

narginchk(4,9);
pdebc = thePde;
if isa(thePde.Mesh, 'pde.FEMesh')
    msh = thePde.Mesh;
    [p,~,t] = msh.meshToPet();
    pdebc = pde.PDEModel(thePde.PDESystemSize);  
    pdebc.BoundaryConditions = thePde.BoundaryConditions;    
    c = varargin{1};
    a = varargin{2};
    f = varargin{3};
    idxstart = 4;
    if ~thePde.IsTwoD
       if nargin == 7 || (nargin == 6 && isscalar(varargin{4}))
         error(message('pde:pdeModel:sdl3D'));
       end
    end
else
    if ~thePde.IsTwoD
      error(message('pde:assema:pdemodelNoMesh'));  
    end
    p = varargin{1};  
    t = varargin{2};  
    c = varargin{3};
    a = varargin{4};
    f = varargin{5};
    idxstart = 6;
end

nrp = size(p,1);
if(nrp==3)
  import pde.internal.pdeEquations;
  if nargin > idxstart
    asm = pdeEquations(pdebc,p,[],t,c,a,f,varargin{idxstart:end});
  else
    asm = pdeEquations(pdebc,p,[],t,c,a,f);  
  end
  try
    [ ks, km, fm ] = getKMF(asm);
  catch ex
    throw(ex);
  end
elseif(nrp==2)
  [ks,km,fm] = pde.internal.assema2DImpl(pdebc, p, t, c, a, f, varargin{idxstart:end});
else
  error(message('pde:assema:PLength', nrp));
end

end

