function nflies_per_roi = CountFliesPerROI(readframe,bgmed,nframes,roidata,roiparams,tracking_params)

nrois = numel(roidata.centerx);

% do background subtraction to count flies in each roi
framessample = round(linspace(1,nframes,roiparams.nframessample));
areassample = cell(nrois,roiparams.nframessample);

for i = 1:roiparams.nframessample,
  im = readframe(framessample(i));
  switch tracking_params.bgmode,
    case 'DARKBKGD',
      dbkgd = imsubtract(im,bgmed);
    case 'LIGHTBKGD',
      dbkgd = imsubtract(bgmed,im);
    case 'OTHERBKGD',
      dbkgd = imabsdiff(im,bgmed);
    otherwise
      error('Unknown background type');
  end
  
  % threshold
  isfore = dbkgd >= tracking_params.bgthresh;
  
  for j = 1:nrois,
    roibb = roidata.roibbs(j,:);
    isforebb = isfore(roibb(3):roibb(4),roibb(1):roibb(2));
    isforebb(~roidata.inrois{j}) = false;
    cc = bwconncomp(isforebb);
    areassample{j,i} = cellfun(@numel,cc.PixelIdxList);    
  end
  
end

% heuristic: if mode ~= 1, use mode
% otherwise use 99th percentile
nflies_per_roi = nan(1,nrois);
for i = 1:nrois,
  nccs = cellfun(@(x) nnz(x >= tracking_params.minccarea),areassample(i,:));
  mode_nccs = mode(nccs);
  if mode_nccs ~= 1,
    nflies_per_roi(i) = mode_nccs;
  else
    nflies_per_roi(i) = round(prctile(nccs,99));
  end
end