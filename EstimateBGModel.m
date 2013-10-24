function [bgmed,bgfixdata] = EstimateBGModel(readframe,nframes,tracking_params,dofixbg)

if nargin < 4,
  dofixbg = false;
end

% compute background model
buffer = readframe(1);
ncolors = size(buffer,3);
if ncolors > 1,
  buffer = rgb2gray(buffer);
end
buffer = repmat(buffer,[1,1,tracking_params.bg_nframes]);
frames = round(linspace(1,min(nframes,tracking_params.bg_lastframe),tracking_params.bg_nframes));
for i = 1:tracking_params.bg_nframes,
  if mod(i,10) == 0,
    fprintf('*');
  end
  t = frames(i);
  im = readframe(t);
  if ncolors > 1,
    im = rgb2gray(im);
  end
  buffer(:,:,i) = im;
end
fprintf('\n');
bgmed = uint8(median(single(buffer),3));

if dofixbg,
  [bgmed,bgfixdata] = FixBgModel(bgmed,moviefile,tracking_params);
else
  bgfixdata = [];
end
