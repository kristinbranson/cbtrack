function roidata = CourtshipBowlDetectROIs(expdir,varargin)

version = '0.1';
timestamp = datestr(now,TimestampFormat);

roidata = struct;
roidata.cbdetectrois_version = version;
roidata.cbdetectrois_timestamp = timestamp;

%% parse inputs
ParseCourtshipBowlParams;
[~,experiment_name] = fileparts(expdir);

% read parameters
params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.roidetectparamsfilestr));
params.roimus = reshape(params.roimus,[2,numel(params.roimus)/2])';
for i = 1:2:numel(leftovers)-1, %#ok<USENS>
  if isfield(params,leftovers{i}),
    params.(leftovers{i}) = leftovers{i+1};
  end
end
tracking_params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.trackingparamsfilestr));

%% open movie
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
imheight = headerinfo.nr;
imwidth = headerinfo.nc;

%% estimate the background model

% compute background model
fprintf('Computing background model for %s...\n',experiment_name);
buffer = readframe(1);
buffer = repmat(buffer,[1,1,tracking_params.bg_nframes]);
frames = round(linspace(1,nframes,tracking_params.bg_nframes));
for i = 1:tracking_params.bg_nframes,
  if mod(i,10) == 0,
    fprintf('*');
  end
  t = frames(i);
  buffer(:,:,i) = readframe(t);
end
fprintf('\n');
bgmed = uint8(median(single(buffer),3));

savefile = fullfile(expdir,dataloc_params.bgmatfilestr);
if exist(savefile,'file'),
  delete(savefile);
end
save(savefile,'bgmed','version','timestamp');

%% detect circles

fprintf('Detecting rois for %s...\n',experiment_name);
[roidata.centerx,roidata.centery,roidata.radii,roidata.scores] = detectCourtshipBowlROIs(bgmed,params);
nrois = numel(roidata.centerx);

%% compute image to real-world transform

roidata.roidiameter_mm = params.roidiameter_mm;
roidata.pxpermm = nanmean(roidata.radii) / (roidata.roidiameter_mm/2);

% find rotations
if isfield(params,'roirows'),
    
  tmp = str2double(params.roirows);
  fnis = find(isnan(tmp));
  params.roirownames = params.roirows(fnis);
  fnis(end+1) = numel(tmp)+1;
  params.roirows = cell(1,numel(fnis)-1);
  for i = 1:numel(fnis)-1,
    params.roirows{i} = tmp(fnis(i)+1:fnis(i+1)-1);
  end
  
  thetas = [];
  for i = 1:numel(params.roirows),
    for roii1 = 1:numel(params.roirows{i})-1,
      roi1 = params.roirows{i}(roii1);
      for roii2 = roii1+1:numel(params.roirows{i}),
        roi2 = params.roirows{i}(roii2);
        thetas(end+1) = modrange(atan2(roidata.centery(roi2)-roidata.centery(roi1),...
          roidata.centerx(roi2)-roidata.centerx(roi1)),-pi/2,pi/2); %#ok<AGROW>
      end
    end
  end
  
  dtheta = modrange(thetas-thetas(1),-pi,pi);
  meantheta = modrange(thetas(1) + mean(dtheta),-pi,pi);
  didaddrotateby = false;
  if isfield(params,'baserotatebyperrigbowl'),
    rigbowl = '';
    metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
    if exist(metadatafile,'file'),
      metadata = ReadMetadataFile(metadatafile);
      if isfield(metadata,'rig') && isfield(metadata,'bowl'),
        rigbowl = sprintf('%d%s',metadata.rig,metadata.bowl);
      end
    end
    if isempty(rigbowl),
      metadata = parseExpDir(expdir,true);
      if isstruct(metadata) && isfield(metadata,'rig') && isfield(metadata,'bowl'),
        rigbowl = sprintf('%d%s',metadata.rig,metadata.bowl);
      end
    end
    if ~isempty(rigbowl),
      i = find(strcmp(params.baserotatebyperrigbowl,rigbowl));
      baserotateby = str2double(params.baserotatebyperrigbowl{i+1});
      meantheta = modrange(-meantheta+baserotateby*pi/180,-pi,pi);
      didaddrotateby = true;
    end
  end
  if ~didaddrotateby && isfield(params,'baserotateby'),
    meantheta = modrange(-meantheta+params.baserotateby*pi/180,-pi,pi);
  end
  roidata.rotateby = -meantheta;
  
end

%% create masks

[XGRID,YGRID] = meshgrid(1:imwidth,1:imheight);
roibbs = [max(1,floor(roidata.centerx(:)-roidata.radii(:))),...
  min(imwidth,ceil(roidata.centerx(:)+roidata.radii(:))),...
  max(1,floor(roidata.centery(:)-roidata.radii(:))),...
  min(imheight,ceil(roidata.centery(:)+roidata.radii(:)))];

idxroi = zeros(imheight,imwidth);
inrois = cell(1,nrois);
for i = 1:nrois,
  bb = roibbs(i,:);
  inrois{i} = (XGRID(bb(3):bb(4),bb(1):bb(2)) - roidata.centerx(i)).^2 + ...
    (YGRID(bb(3):bb(4),bb(1):bb(2)) - roidata.centery(i)).^2 ...
    <= roidata.radii(i)^2;
  tmp = idxroi(bb(3):bb(4),bb(1):bb(2));
  tmp(inrois{i}) = i;
  idxroi(bb(3):bb(4),bb(1):bb(2)) = tmp;
end
roidata.roibbs = roibbs;
roidata.idxroi = idxroi;
roidata.inrois = inrois;

%% count number of flies in each roi

roidata.nflies_per_roi = CountFliesPerROI(readframe,bgmed,nframes,roidata,params,tracking_params);

%% save results

savefile = fullfile(expdir,dataloc_params.roidatamatfilestr);
if exist(savefile,'file'),
  delete(savefile);
end
save(savefile,'-struct','roidata');
  
%% plot results

hfig = 18159;
figure(hfig);
clf;
figpos = [10,10,imwidth,imheight];
set(hfig,'Units','pixel','Position',figpos);
axes('Position',[0,0,1,1]);
imagesc(readframe(round(nframes/2)),[0,255]);
axis image;
axis off;
colormap gray;
hold on;
colors = jet(nrois)*.7;
for i = 1:nrois,
  drawellipse(roidata.centerx(i),roidata.centery(i),0,roidata.radii(i),roidata.radii(i),'Color',colors(i,:));
  text(roidata.centerx(i),roidata.centery(i),{sprintf('%d flies',roidata.nflies_per_roi(i)),sprintf('score = %d',round(roidata.scores(i)))},...
    'Color',colors(i,:),'HorizontalAlignment','center','VerticalAlignment','middle');
end

set(hfig,'Units','pixels','Position',figpos);
imsavename = fullfile(expdir,dataloc_params.roiimagefilestr);
if exist(imsavename,'file'),
  delete(imsavename);
end
save2png(imsavename,hfig);

%% clean up

% close movie
if fid > 0,
  try
    fclose(fid);
  catch ME,
    warning('Could not close movie file: %s',getReport(ME));
  end
end

if isdeployed && ishandle(hfig),
  delete(hfig);
end