function gm = analyticToDiscrete(self)
% analyticToDiscrete - Create discrete geometry from 2D analytic geometry
%
% Internal use only

% Copyright 2017 The MathWorks, Inc.
    gm = [];
    if isempty(self.Geometry)
        gm = [];
        return;
    end    
    if self.IsTwoD && isa(self.Geometry, 'pde.AnalyticGeometry')         
        F = self.facetAnalyticGeometry(self.Geometry);
        gm = analyticToDiscrete(F);
    end
end