function roidata = CourtshipBowlDetectROIs(expdir,varargin)

version = '0.1.3';
timestamp = datestr(now,TimestampFormat);

roidata = struct;
roidata.cbdetectrois_version = version;
roidata.cbdetectrois_timestamp = timestamp;

%% parse inputs
ParseCourtshipBowlParams;
[dofixbg,leftovers] = myparse_nocheck(leftovers,'dofixbg',false); %#ok<NODEF>
[~,experiment_name] = fileparts(expdir);

params = cbparams.detect_rois;
tracking_params = cbparams.track;

% reformat roimus
params.roimus = [params.roimus.x(:),params.roimus.y(:)];

% read parameters
for i = 1:2:numel(leftovers)-1, 
  if isfield(params,leftovers{i}),
    params.(leftovers{i}) = leftovers{i+1};
  end
end
roidata.params = cbparams;

%% open log file

if isfield(cbparams.dataloc,'roi_log'),
  logfile = fullfile(expdir,cbparams.dataloc.roi_log.filestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

fprintf(logfid,'\n\n***\nRunning CourtshipBowlDetectROIs version %s analysis_protocol %s at %s\n',version,real_analysis_protocol,timestamp);
roidata.analysis_protocol = real_analysis_protocol;

%% open movie

moviefile = fullfile(expdir,cbparams.dataloc.movie.filestr);
fprintf(logfid,'Opening movie file %s...\n',moviefile);
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
im = readframe(1);
[imheight,imwidth,~] = size(im);

%% estimate the background model

% compute background model
fprintf(logfid,'Computing background model for %s...\n',experiment_name);
buffer = readframe(1);
buffer = repmat(buffer,[1,1,tracking_params.bg_nframes]);
frames = round(linspace(1,min(nframes,tracking_params.bg_lastframe),tracking_params.bg_nframes));
for i = 1:tracking_params.bg_nframes,
  if mod(i,10) == 0,
    fprintf('*');
  end
  t = frames(i);
  buffer(:,:,i) = readframe(t);
end
fprintf('\n');
bgmed = uint8(median(single(buffer),3));

if dofixbg,
  [bgmed,bgfixdata] = FixBgModel(bgmed,moviefile,tracking_params);
end

savefile = fullfile(expdir,cbparams.dataloc.bgmat.filestr);
fprintf(logfid,'Saving background model to file %s...\n',savefile);
if exist(savefile,'file'),
  delete(savefile);
end
save(savefile,'bgmed','version','timestamp','tracking_params');

bgimagefile = fullfile(expdir,cbparams.dataloc.bgimage.filestr);
fprintf(logfid,'Saving image of background model to file %s...\n',bgimagefile);
imwrite(bgmed,bgimagefile,'png');

%% detect circles

fprintf(logfid,'Detecting rois for %s...\n',experiment_name);
[roidata.centerx,roidata.centery,roidata.radii,roidata.scores] = DetectCourtshipBowlROIs(bgmed,params);
nrois = numel(roidata.centerx);
fprintf(logfid,'Detected %d ROIs with mean radius %f (std = %f, min = %f, max = %f)\n',...
  nrois,mean(roidata.radii),std(roidata.radii,1),min(roidata.radii),max(roidata.radii));

%% compute image to real-world transform

roidata.roidiameter_mm = params.roidiameter_mm;
roidata.pxpermm = nanmean(roidata.radii) / (roidata.roidiameter_mm/2);
fprintf(logfid,'Computed pxpermm = %f\n',roidata.pxpermm);

% find rotations
if isfield(params,'roirows'),
    
  roirownames = fieldnames(params.roirows);
  
  thetas = [];
  for i = 1:numel(roirownames),
    roiscurr = params.roirows.(roirownames{i});
    for roii1 = 1:numel(roiscurr)-1,
      roi1 = roiscurr(roii1);
      for roii2 = roii1+1:numel(roiscurr),
        roi2 = roiscurr(roii2);
        thetas(end+1) = modrange(atan2(roidata.centery(roi2)-roidata.centery(roi1),...
          roidata.centerx(roi2)-roidata.centerx(roi1)),-pi/2,pi/2); %#ok<AGROW>
      end
    end
  end
  
  dtheta = modrange(thetas-thetas(1),-pi,pi);
  meantheta = modrange(thetas(1) + mean(dtheta),-pi,pi);
  
  fprintf(logfid,'Based on ROI centroids, we want to rotate by %f deg (std = %f, min = %f, max = %f)\n',...
    meantheta*180/pi,std(dtheta,1)*180/pi,modrange(thetas(1) + min(dtheta),-pi,pi)*180/pi,...
    modrange(thetas(1) + max(dtheta),-pi,pi)*180/pi);
  
  didaddrotateby = false;
  if isfield(params,'baserotatebyperrigbowl'),
    rigbowl = '';
    metadatafile = fullfile(expdir,cbparams.dataloc.metadata.filestr);
    if exist(metadatafile,'file'),
      metadata = ReadMetadataFile(metadatafile);
      if isfield(metadata,'rig') && isfield(metadata,'bowl'),
        rigbowl = sprintf('rig%dbowl%s',metadata.rig,metadata.bowl);
      end
    end
    if isempty(rigbowl),
      metadata = parseExpDir(expdir,true);
      if isstruct(metadata) && isfield(metadata,'rig') && isfield(metadata,'bowl'),
        rigbowl = sprintf('rig%dbowl%s',metadata.rig,metadata.bowl);
      end
    end
    if ~isempty(rigbowl) && isfield(params.baserotatebyperrigbowl,rigbowl),
      baserotateby = params.baserotatebyperrigbowl.(rigbowl);
      meantheta = modrange(-meantheta+baserotateby*pi/180,-pi,pi);
      didaddrotateby = true;
      fprintf(logfid,'Adding in baserotateby %f selected for %s\n',baserotateby,rigbowl);
    end
  end
  if ~didaddrotateby && isfield(params,'baserotateby'),
    meantheta = modrange(-meantheta+params.baserotateby*pi/180,-pi,pi);
    fprintf(logfid,'Adding in default baserotateby %f\n',params.baserotateby);
  end
  roidata.rotateby = -meantheta;
  fprintf(logfid,'Final rotateby = %f\n',roidata.rotateby*180/pi);
  
end

%% create masks

fprintf(logfid,'Creating ROI masks...\n');

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

fprintf(logfid,'Counting number of flies in each ROI...\n');

roidata.nflies_per_roi = CountFliesPerROI(readframe,bgmed,nframes,roidata,rigbowl,params,tracking_params);

fprintf(logfid,'nflies\tnrois\n');
for i = 0:2,
  fprintf(logfid,'%d\t%d\n',i,nnz(roidata.nflies_per_roi==i));
end
fprintf(logfid,'>2\t%d\n',nnz(roidata.nflies_per_roi>2));
fprintf(logfid,'ignored\t%d\n',nnz(isnan(roidata.nflies_per_roi)));
fprintf(logfid,'\n');


%% save results

savefile = fullfile(expdir,cbparams.dataloc.roidatamat.filestr);
fprintf(logfid,'Saving results to file %s...\n',savefile);
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
  if isnan(roidata.nflies_per_roi(i)),
    text(roidata.centerx(i),roidata.centery(i),{'ignored',sprintf('score = %d',round(roidata.scores(i)))},...
      'Color',colors(i,:),'HorizontalAlignment','center','VerticalAlignment','middle');
  else
    text(roidata.centerx(i),roidata.centery(i),{sprintf('%d flies',roidata.nflies_per_roi(i)),sprintf('score = %d',round(roidata.scores(i)))},...
      'Color',colors(i,:),'HorizontalAlignment','center','VerticalAlignment','middle');
  end
end

set(hfig,'Units','pixels','Position',figpos);
imsavename = fullfile(expdir,cbparams.dataloc.roiimage.filestr);
fprintf(logfid,'Outputting visualization fo results to %s...\n',imsavename);
if exist(imsavename,'file'),
  delete(imsavename);
end
save2png(imsavename,hfig);

%% clean up

fprintf(logfid,'Cleaning up...\n');

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

fprintf(logfid,'Finished running CourtshipBowlDetectROIS at %s.\n',datestr(now,'yyyymmddTHHMMSS'));
if logfid > 1,
  fclose(logfid);
end