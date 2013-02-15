function [trxcurr,pred] = TrackTwoFliesOneFrame(dbkgd,isfore,pred,trxprev,roidata,params)

trxcurr = trxprev;
% loop over rois  
for roii = 1:roidata.nrois,
  
  if roidata.nflies_per_roi(roii) == 0,
    continue;
  end
  
  % crop out roi  
  roibb = roidata.roibbs(roii,:);
  isforebb = isfore(roibb(3):roibb(4),roibb(1):roibb(2));
  dbkgdbb = double(dbkgd(roibb(3):roibb(4),roibb(1):roibb(2)));
  dbkgdbb(~roidata.inrois{roii}) = 0;
  isforebb(~roidata.inrois{roii}) = false;
  
  % main work
  [trxcurr(roii),pred(roii)] = ...
    TrackTwoFliesOneFrameOneROI(isforebb,dbkgdbb,pred(roii),trxcurr(roii),roidata.nflies_per_roi(roii),params);

end