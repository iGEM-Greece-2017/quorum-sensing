function [x,y]= bactDomainGeom(x0,bactSize,domainLim,args)
% Defines a geometry consisting of a hexagonal domain with many elliptical bacteria inside.
% Input:
% - x0 {nBact,2}: Center of each bacterium
% - bactSize {2}: Bacteria size (same for all)
% - domainLim: Radius of the hexagonal domain
% - args: cell array of the extra arguments needed by geometry functions

nBact= size(x0,1);
switch length(args)
  case 0
    x= 6+2*nBact;           % Line segments
    
  case 1
    bs= args{1};
    A= zeros(4,length(bs));
    % Domain
    edges= bs<=6;
    if any(edges(:)), A(:,edges)= util.hexFun(domainLim, {bs(edges)}); end
    % Bacteria
    for b= 1:nBact
      lims= [6+2*(b-1), 6+2*b];
      edges= bs>lims(1) & bs<=lims(2);
      if any(edges(:))
        bsLocal= bs(edges) -6-(b-1)*2;
        A(:,edges)= util.ellipseFun(x0(b,:), bactSize, {bsLocal});
        % Globalize parameter values.
        % TODO: The transformation (amount added) should be parameterized by the final param
        %   value returned by each geometry
        A(1:2,edges)= A(1:2,edges) +2+(b-1)*2*pi;
        A(3,edges)= b+1;  % Assign face numbers
        A(4,edges)= 1;
      end
    end
    x= A;
    
  case 2
    bs= args{1};
    t= args{2};
    x= zeros(size(t)); y= zeros(size(t));
    if numel(bs) > 1, assert(all(size(bs) == size(t))); end
    % Domain
    edges= bs<=6;
    if any(edges(:))
      tsel= edges;
      if numel(tsel)==1, tsel= 1:length(t); end
      [x(tsel),y(tsel)]= util.hexFun(domainLim, {bs(edges), t(tsel)});
    end
    % Bacteria
    for b= 1:nBact
      lims= [6+(b-1)*2, 6+b*2];
      edges= bs>lims(1) & bs<=lims(2);
      if any(edges(:))
        tsel= edges;
        if numel(tsel)==1, tsel= 1:length(t); end
        bsLocal= bs(edges) -6-(b-1)*2;
        tLocal= t(tsel)    -2-(b-1)*2*pi;
        [x(tsel),y(tsel)]= util.ellipseFun(x0(b,:), bactSize, {bsLocal, tLocal});
      end
    end
end
