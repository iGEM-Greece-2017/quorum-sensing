function [x,y]= bactAgarGeom(x0,bactSize,domainLim,args)
% Define a rectangular domain, horizontal, with its upper side at y=0 and
% extending downwards, with bacteria embedded inside. Represents a vertical cut
% of a cylindrical agar plate: y->z and x->r. Implies axial symmetry.

nBact= size(x0,1);
switch length(args)
  case 0
    x= 4+2*nBact;
  case 1
    bs= args{1};
    A= zeros(4,length(bs));
    l= [2,1].*domainLim;
    % Domain
    edges= bs<=4;
    if any(edges(:)), A(:,edges)= util.rectangleFun([-domainLim(1),0],l, {bs(edges)}); end
    % Bacteria
    for b= 1:nBact
      lims= [4+2*(b-1), 4+2*b];
      edges= bs>lims(1) & bs<=lims(2);
      if any(edges(:))
        bsLocal= bs(edges) -4-(b-1)*2;
        A(:,edges)= util.ellipseFun(x0(b,:), bactSize, {bsLocal});
        % Globalize parameter values.
        % TODO: The transformation (amount added) should be parameterized by the final param
        %   value returned by each geometry
        A(1:2,edges)= A(1:2,edges) +2*sum(l)+(b-1)*2*pi;
        A(3,edges)= b+1;  % Assign face numbers
        A(4,edges)= 1;
      end
    end
    x= A;
    
  case 2
    bs= args{1};
    t= args{2};
    l= [2,1].*domainLim;
    x= zeros(size(t)); y= zeros(size(t));
    if numel(bs) > 1, assert(all(size(bs) == size(t))); end
    % Domain
    edges= bs<=4;
    if any(edges(:))
      tsel= edges;
      if numel(tsel)==1, tsel= 1:length(t); end
      [x(tsel),y(tsel)]= util.rectangleFun([-domainLim(1),0],l, {bs(edges), t(tsel)});
    end
    % Bacteria
    for b= 1:nBact
      lims= [4+(b-1)*2, 4+b*2];
      edges= bs>lims(1) & bs<=lims(2);
      if any(edges(:))
        tsel= edges;
        if numel(tsel)==1, tsel= 1:length(t); end
        bsLocal= bs(edges) -4-(b-1)*2;
        tLocal= t(tsel)    -2*sum(l)-(b-1)*2*pi;
        [x(tsel),y(tsel)]= util.ellipseFun(x0(b,:), bactSize, {bsLocal, tLocal});
      end
    end
end
