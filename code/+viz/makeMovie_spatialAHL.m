function makeMovie_spatialAHL(hist,initTime, filename)

makeVideo= nargin==3;
totalAHLlims= [arrayfun(@(x) min(x.AHL(:)), hist), arrayfun(@(x) max(x.AHL(:)), hist)];
totalAHLlims= [min(totalAHLlims(:,1)), max(totalAHLlims(:,2))];
steps= length(hist);
ahl_g= imagesc(hist(1).AHL, totalAHLlims);
title(num2str(initTime));

% Generate movie frames
tic;
if makeVideo, frames(steps)= struct('cdata',[],'colormap',[]); end;
for i=1:steps
  ahl_g.CData= hist(i).AHL;
  ahl_g.Parent.Title.Text.String= num2str(i+initTime-1);
  drawnow;
  if makeVideo, frames(i)= getframe; end;
end
toc
% Write movie to file
if makeVideo
  tic;
  vw= VideoWriter(filename, 'Archival');
  vw.FrameRate= 25;
  vw.open();
  vw.writeVideo(frames(1:2:end));
  vw.close();
  toc
end
