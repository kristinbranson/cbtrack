function trackdata = TrackTwoFlies(moviefile,bgmed,roidata,params,varargin)

version = '0.1.1';
timestamp = datestr(now,TimestampFormat);
trackdata = struct;
trackdata.tracktwoflies_version = version;
trackdata.tracktwoflies_timestamp = timestamp;

tmpfilename = sprintf('TmpResultsTrackTwoFlies_%s.mat',datestr(now,TimestampFormat));
[restart,tmpfilename] = myparse(varargin,'restart','','tmpfilename',tmpfilename);

dorestart = false;
if ~isempty(restart),
  fprintf('Loading temporary results from file %s...\n',restart);
  load(restart);
  restartstage = stage; %#ok<NODEF>
  dorestart = true;
end
fprintf('TrackTwoFlies temporary results saved to file %s\n',tmpfilename);

SetBackgroundTypes;
flycolors = {'r','b'};
stages = {'maintracking','chooseorientations','trackwings1','assignids','chooseorientations2','trackwings2'};

%% open movie

fprintf('Opening movie...\n');
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);

%% initialize

nrois = roidata.nrois;
nframes_track = min(params.lastframetrack,nframes)-params.firstframetrack+1;

if ~dorestart,
  fprintf('Allocating...\n');
  trxx = nan(2,nrois,nframes_track);
  trxy = nan(2,nrois,nframes_track);
  trxa = nan(2,nrois,nframes_track);
  trxb = nan(2,nrois,nframes_track);
  trxtheta = nan(2,nrois,nframes_track);
  trxarea = nan(2,nrois,nframes_track);
  istouching = nan(nrois,nframes_track);
  gmm_isbadprior = nan(nrois,nframes_track);
  
  pred = struct;
  pred.mix = gmm(2,2,'full');
  pred.x = nan(2,1);
  pred.y = nan(2,1);
  pred.theta = nan(2,1);
  pred.area = nan(2,1);
  pred.isfirstframe = true;
  pred = repmat(pred,[1,nrois]);
  
  trxcurr = struct;
  trxcurr.x = nan(2,1);
  trxcurr.y = nan(2,1);
  trxcurr.a = nan(2,1);
  trxcurr.b = nan(2,1);
  trxcurr.theta = nan(2,1);
  trxcurr.area = nan(2,1);
  trxcurr.istouching = nan;
  trxcurr.gmm_isbadprior = nan;
  trxcurr = repmat(trxcurr,[1,roidata.nrois]);
  
end

hell = nan(2,nrois);
htrx = nan(2,nrois);
him = nan;

%% loop over frames

stage = 'maintracking'; 

if ~dorestart || find(strcmp(stage,stages)) >= find(strcmp(restartstage,stages)),

if dorestart && strcmp(stage,restartstage),
  startframe = t; %#ok<NODEF>
else
  startframe = params.firstframetrack;
end

fprintf('Main tracking...\n');
for t = startframe:min(params.lastframetrack,nframes),
  
  iframe = t - params.firstframetrack + 1;
  
  if mod(iframe,1000) == 0,
    fprintf('Frame %d / %d\n',iframe,nframes_track);
  end
  
  % read in frame
  im = readframe(t);
  
  % subtract off background
  switch params.bgmode,
    case DARKBKGD,
      dbkgd = imsubtract(im,bgmed);
    case LIGHTBKGD,
      dbkgd = imsubtract(bgmed,im);
    case OTHERBKGD,
      dbkgd = imabsdiff(im,bgmed);
  end

  % threshold
  isfore = dbkgd >= params.bgthresh;

  % track this frame
  trxprev = trxcurr;
  [trxcurr,pred] = TrackTwoFliesOneFrame(dbkgd,isfore,pred,trxprev,roidata,params);

  trxx(:,:,iframe) = cat(2,trxcurr.x);
  trxy(:,:,iframe) = cat(2,trxcurr.y);
  trxa(:,:,iframe) = cat(2,trxcurr.a);
  trxb(:,:,iframe) = cat(2,trxcurr.b);
  trxtheta(:,:,iframe) = cat(2,trxcurr.theta);
  trxarea(:,:,iframe) = cat(2,trxcurr.area);
  istouching(:,iframe) = cat(1,trxcurr.istouching);
  gmm_isbadprior(:,iframe) = cat(1,trxcurr.gmm_isbadprior);
  
  % plot
  if params.DEBUG,
    
    if t == params.firstframetrack || ~ishandle(him),
      hold off;
      him = imagesc(im,[0,255]);
      axis image;
      colormap gray;
      hold on;
      isnewplot = true;
    else
      set(him,'CData',im);
      isnewplot = false;
    end
    title(num2str(t));
    
    for roii = 1:roidata.nrois,
      
      roibb = roidata.roibbs(roii,:);
      offx = roibb(1)-1;
      offy = roibb(3)-1;
      
      for i = 1:2,
        if isnewplot || ~ishandle(hell(i,roii)),
          hell(i,roii) = drawellipse(trxx(i,roii,iframe)+offx,trxy(i,roii,iframe)+offy,trxtheta(i,roii,iframe),trxa(i,roii,iframe),trxb(i,roii,iframe),[flycolors{i},'-']);
        else
          updateellipse(hell(i,roii),trxx(i,roii,iframe)+offx,trxy(i,roii,iframe)+offy,trxtheta(i,roii,iframe),trxa(i,roii,iframe),trxb(i,roii,iframe));
        end
        if isnewplot || ~ishandle(htrx(i,roii)),
          htrx(i,roii) = plot(squeeze(trxx(i,roii,max(iframe-30,1):iframe)+offx),squeeze(trxy(i,roii,max(iframe-30,1):iframe)+offy),[flycolors{i},'.-']);
        else
          set(htrx(i,roii),'XData',squeeze(trxx(i,roii,max(iframe-30,1):iframe)+offx),...
            'YData',squeeze(trxy(i,roii,max(iframe-30,1):iframe)+offy));
        end
      end
    end
    
  end
  
  if params.DEBUG || mod(t,10) == 0,
    drawnow;
  end
  
  if mod(iframe,5000) == 0,
    save(tmpfilename,'trxx','trxy','trxa','trxb','trxtheta','trxarea','istouching','gmm_isbadprior','pred','trxcurr','t','params','moviefile','bgmed','roidata','stage');
  end
  
end

save(tmpfilename,'trxx','trxy','trxa','trxb','trxtheta','trxarea','istouching','gmm_isbadprior','pred','trxcurr','t','params','moviefile','bgmed','roidata','stage');

%% correct for bounding box of rois

fprintf('Correcting for ROI bounding boxes...\n');
for roii = 1:nrois,
  roibb = roidata.roibbs(roii,:);
  trxx(:,roii,:) = trxx(:,roii,:) + roibb(1) - 1;
  trxy(:,roii,:) = trxy(:,roii,:) + roibb(3) - 1;
end

%% reformat

fprintf('Reformatting...\n');
if isfield(headerinfo,'timestamps'),
  timestamps = headerinfo.timestamps;
elseif isfield(headerinfo,'FrameRate'),
  timestamps = (0:nframes-1)/headerinfo.FrameRate;
elseif isfield(headerinfo,'fps'),
  timestamps = (0:nframes-1)/headerinfo.fps;
else
  warning('No frame rate info found for movie');
  timestamps = nan(1,nframes);
end
trackdata.timestamps = timestamps;

trackdata.trx = struct;
j = 1;
fly2roiid = [];
for i = 1:nrois,
  for jj = 1:roidata.nflies_per_roi(i),
    trackdata.trx(j).x = reshape(trxx(jj,i,:),[1,nframes_track]);
    trackdata.trx(j).y = reshape(trxy(jj,i,:),[1,nframes_track]);
    % saved quarter-major, quarter-minor axis
    trackdata.trx(j).a = reshape(trxa(jj,i,:)/2,[1,nframes_track]);
    trackdata.trx(j).b = reshape(trxb(jj,i,:)/2,[1,nframes_track]);
    trackdata.trx(j).theta = reshape(trxtheta(jj,i,:),[1,nframes_track]);

    trackdata.trx(j).firstframe = params.firstframetrack;
    trackdata.trx(j).endframe = params.firstframetrack+nframes_track-1;
    trackdata.trx(j).nframes = nframes_track;
    trackdata.trx(j).off = 1-params.firstframetrack;
    trackdata.trx(j).roi = i;
    trackdata.trx(j).arena = struct;
    trackdata.trx(j).arena.arena_radius_mm = roidata.radii(i);
    trackdata.trx(i).arena.arena_center_mm_x = roidata.centerx(i);
    trackdata.trx(i).arena.arena_center_mm_y = roidata.centery(i);
    %trackdata.trx(j).roipts = rois{i};
    %trackdata.trx(j).roibb = roibbs(i,:);
    trackdata.trx(j).moviefile = moviefile;
    trackdata.trx(j).dt = diff(timestamps(trackdata.trx(j).firstframe:trackdata.trx(j).endframe));
    trackdata.trx(j).timestamps = timestamps(trackdata.trx(j).firstframe:trackdata.trx(j).endframe);    
    fly2roiid(j) = jj; %#ok<AGROW>
    j = j + 1;
  end
end
trackdata.istouching = istouching;
trackdata.gmm_isbadprior = gmm_isbadprior;

nflies = numel(trackdata.trx);

end

%% resolve head/tail ambiguity

fprintf('Choosing orientations 1...\n');
stage = 'chooseorientations'; 
save(tmpfilename,'trackdata','params','moviefile','bgmed','roidata','nflies','fly2roiid','stage');

if ~dorestart || find(strcmp(stage,stages)) >= find(strcmp(restartstage,stages)),

for i = 1:nflies,
  x = trackdata.trx(i).x;
  y = trackdata.trx(i).y;
  theta = trackdata.trx(i).theta;
  trackdata.trx(i).theta = choose_orientations(x,y,theta,params.choose_orientations_velocity_angle_weight,params.choose_orientations_max_velocity_angle_weight);
end

end
  
%% plot updated results

if params.DEBUG > 1,

  hell = nan(1,nflies);
  htrx = nan(1,nflies);
  him = nan;

  for iframe = 1:min(nframes_track,500),
    
    t = iframe+params.firstframetrack-1;
    im = readframe(t);
    
    if t == params.firstframetrack || ~exist('him','var') || ~ishandle(him),
      hold off;
      him = imagesc(im,[0,255]);
      axis image;
      colormap gray;
      hold on;
      isnewplot = true;
    else
      set(him,'CData',im);
      isnewplot = false;
    end
    title(num2str(t));

    for i = 1:nflies,
      fly = fly2roiid(i);
      if isnewplot || ~ishandle(hell(i)),
        hell(i) = drawflyo(trackdata.trx(i).x(iframe),trackdata.trx(i).y(iframe),...
          trackdata.trx(i).theta(iframe),trackdata.trx(i).a(iframe),...
          trackdata.trx(i).b(iframe),[flycolors{fly},'-']);
      else
        updatefly(hell(i),trackdata.trx(i).x(iframe),trackdata.trx(i).y(iframe),trackdata.trx(i).theta(iframe),...
          trackdata.trx(i).a(iframe),trackdata.trx(i).b(iframe));
      end
      if isnewplot || ~ishandle(htrx(i)),
        htrx(i) = plot(trackdata.trx(i).x(max(iframe-30,1):iframe),trackdata.trx(i).y(max(iframe-30,1):iframe),[flycolors{fly},'.-']);
      else
        set(htrx(i),'XData',trackdata.trx(i).x(max(iframe-30,1):iframe),...
          'YData',trackdata.trx(i).y(max(iframe-30,1):iframe));
      end
    end
    drawnow;
  end
  
end

%% track wings if using wings

stage = 'trackwings1'; 
save(tmpfilename,'trackdata','params','moviefile','bgmed','roidata','nflies','fly2roiid','stage');

didtrackwings = false;

if ~dorestart || find(strcmp(stage,stages)) >= find(strcmp(restartstage,stages)),

if strcmp(params.assignidsby,'wingsize'),

  fprintf('Tracking wings 1...\n');
  
  [nr,nc,~] = size(readframe(1));
  isarena = false(nr,nc);
  [XGRID,YGRID] = meshgrid(1:nc,1:nr);
  for roii = 1:nrois,
    if roidata.nflies_per_roi(roii) == 0,
      continue;
    end
    isarena = isarena | ...
      ( ((XGRID - roidata.centerx(roii)).^2 + ...
      (YGRID - roidata.centery(roii)).^2) ...
      <= roidata.radii(roii)^2 );
  end
  
  [wingtrx,wingperframedata,wingtrackinfo,wingperframeunits] = TrackWingsHelper(trackdata.trx,moviefile,double(bgmed),isarena,params.wingtracking_params,...
    'firstframe',params.firstframetrack,...
    'debug',params.DEBUG);
  didtrackwings = true;
  trackdata.trackwings_timestamp = wingtrackinfo.trackwings_timestamp;
  trackdata.trackwings_version = wingtrackinfo.trackwings_version;
  trackdata.trx = wingtrx;
  trackdata.perframedata = wingperframedata;
  trackdata.perframeunits = wingperframeunits;

end

end

%% assign identities based on size of something

fprintf('Assigning identities based on sizes...\n');

stage = 'assignids'; 
save(tmpfilename,'trackdata','params','moviefile','bgmed','roidata','nflies','fly2roiid','didtrackwings','stage');

if ~dorestart || find(strcmp(stage,stages)) >= find(strcmp(restartstage,stages)),

assignids_nflips = nan(1,nrois);
switch params.assignidsby,
  case 'size',
    mudatafit = nan(2,3,nrois);
    sigmadatafit = nan(2,3,nrois);
  case 'wingsize',
    mudatafit = nan(2,1,nrois);
    sigmadatafit = nan(2,1,nrois);
    wingarea = cell(1,nflies);
    for i = 1:nflies,
      wingarea{i} = trackdata.perframedata.wing_areal{i}+trackdata.perframedata.wing_arear{i};
    end
  otherwise,
    error('Unknown assignidsby value');
end
niters_assignids_em = nan(1,nrois);
cost_assignids = nan(1,nrois);
sigmamotionfit = nan(1,nrois);
idsfit = nan(2,nrois,nframes_track);

oldtrx = trackdata.trx;
if isfield(trackdata,'perframedata'),
  perframefns = fieldnames(trackdata.perframedata);
  oldperframedata = trackdata.perframedata;
end
trxfns = intersect({'x','y','a','b','theta','xwingl','ywingl','xwingr','ywingr','wing_anglel','wing_angler'},fieldnames(trackdata.trx));
for roii = 1:nrois,
  
  if roidata.nflies_per_roi(roii) < 2,
    continue;
  end
  
  flies = find([trackdata.trx.roi]==roii);
  
  x = cat(1,trackdata.trx(flies).x);
  y = cat(1,trackdata.trx(flies).y);
  appearanceweight = double(~trackdata.istouching(roii,:));
  
  switch params.assignidsby,
    case 'size',
      a = cat(1,trackdata.trx(flies).a);
      b = cat(1,trackdata.trx(flies).b);
      area = a.*b.*pi*4;
      iddata = cat(3,area,a,b);
    case 'wingsize',
      iddata = cat(1,wingarea{flies});
    otherwise,
      error('Unknown assignidsby value');
  end
  [idsfit_curr,mudatafit_curr,sigmadatafit_curr,sigmamotionfit_curr,cost_curr,niters_curr] = ...
    AssignIdentities(x,y,iddata,'vel_dampen',params.err_dampen_pos,'appearanceweight',appearanceweight);
  assignids_nflips(roii) = nnz(idsfit_curr(1,1:end-1)~=idsfit_curr(1,2:end));
  fprintf('Roi %d, flipped ids %d times\n',roii,assignids_nflips(roii));
  
  for i = 1:2,
    for j = 1:2,
      idx = idsfit_curr(i,:)==j;
      for k = 1:numel(trxfns),
        trackdata.trx(flies(i)).(trxfns{k})(idx) = oldtrx(flies(j)).(trxfns{k})(idx);
      end
      if isfield(trackdata,'perframedata'),
        for k = 1:numel(perframefns),
          trackdata.perframedata.(perframefns{k}){flies(i)}(idx) = ...
            oldperframedata.(perframefns{k}){flies(j)}(idx);
        end
      end
    end
  end

  mudatafit(:,:,roii) = mudatafit_curr;
  sigmadatafit(:,:,roii) = sigmadatafit_curr;
  niters_assignids_em(roii) = niters_curr;
  cost_assignids(roii) = cost_curr;
  sigmamotionfit(roii) = sigmamotionfit_curr;  
  idsfit(:,roii,:) = idsfit_curr;
end

idx = find(roidata.nflies_per_roi == 2);
areas_male = mudatafit(1,1,idx);
areas_female = mudatafit(2,1,idx);
area_thresh = (mean(areas_male)+mean(areas_female))/2;

idx = find(roidata.nflies_per_roi == 1);
typeperroi = cell(1,nrois);
for ii = 1:numel(idx),
  i = idx(ii);
  fly = find([trackdata.trx.roi]==roii);
  switch params.assignidsby,
    case 'size',
      a = cat(1,trackdata.trx(fly).a);
      b = cat(1,trackdata.trx(fly).b);
      areacurr = a.*b.*pi*4;
    case 'wingsize',
      areacurr = wingarea{fly};
    otherwise,
      error('Unknown assignidsby value');
  end
  meanareacurr = nanmean(areacurr(:));
  fprintf('%d: %f\n',i,meanareacurr);
  if meanareacurr <= area_thresh, 
    typeperroi{i} = params.typesmallval;
  else
    typeperroi{i} = params.typebigval;
  end
end

trackdata.assignids = struct;
trackdata.assignids.nflips = assignids_nflips;
trackdata.assignids.mudatafit = mudatafit;
trackdata.assignids.sigmadatafit = sigmadatafit;
trackdata.assignids.niters_em = niters_assignids_em;
trackdata.assignids.cost = cost_assignids;
trackdata.assignids.sigmamotionfit = sigmamotionfit;
trackdata.assignids.idsfit = idsfit;

for i = 1:nflies,
  roii = trackdata.trx(i).roi;
  if roidata.nflies_per_roi(roii) == 2,
    if fly2roiid(i) == 1,
      trackdata.trx(i).(params.typefield) = repmat({params.typesmallval},[1,nframes_track]);
    else
      trackdata.trx(i).(params.typefield) = repmat({params.typebigval},[1,nframes_track]);
    end
  elseif roidata.nflies_per_roi(i) == 1,
    trackdata.trx(i).(params.typefield) = repmat(typeperroi(i),[1,nframes_track]);
  end
end

end

%% resolve head/tail ambiguity again, since it is pretty quick

stage = 'chooseorientations2'; 
save(tmpfilename,'trackdata','params','moviefile','bgmed','roidata','nflies','fly2roiid','didtrackwings','stage');

if ~dorestart || find(strcmp(stage,stages)) >= find(strcmp(restartstage,stages)),

fprintf('Choosing orientations 2...\n');
isflip = false(nflies,nframes_track);

for i = 1:nflies,
  
  % if there is some kind of flip
  roii = trackdata.trx(i).roi;
  if trackdata.assignids.nflips(roii) == 0,
    continue;
  end
  x = trackdata.trx(i).x;
  y = trackdata.trx(i).y;
  theta = trackdata.trx(i).theta;
  fprintf('Re-choosing orientations for fly %d (nidflips = %d)\n',i,trackdata.assignids.nflips(roii));
  trackdata.trx(i).theta = choose_orientations(x,y,theta,params.choose_orientations_velocity_angle_weight,params.choose_orientations_max_velocity_angle_weight);
  isflip(i,:) = round(abs(modrange(theta-trackdata.trx(i).theta,-pi,pi))/pi) > 0;
  fprintf('N. orientation flips = %d\n',nnz(isflip(i,:)));
  
end

end

%% track wings

stage = 'trackwings2'; 
save(tmpfilename,'trackdata','params','moviefile','bgmed','roidata','nflies','fly2roiid','didtrackwings','isflip','stage');

if ~dorestart || find(strcmp(stage,stages)) >= find(strcmp(restartstage,stages)),

if didtrackwings,
  framestrack = cell(1,nrois);
  for roii = 1:nrois,
    flies = find([trackdata.trx.roi]==roii);
    if isempty(flies),
      continue;
    end
    framestrack{roii} = find(any(isflip(flies,:),1));
  end
  roistrack = find(~cellfun(@isempty,framestrack));
else
  roistrack = 1:nrois;
end
  
if params.dotrackwings && ~isempty(roistrack),
  
  fprintf('Tracking wings 2...\n');
  
  [nr,nc,~] = size(readframe(1));
  isarena = false(nr,nc);
  [XGRID,YGRID] = meshgrid(1:nc,1:nr);
  for roii = roistrack,
    if roidata.nflies_per_roi(roii) == 0,
      continue;
    end
    isarena = isarena | ...
      ( ((XGRID - roidata.centerx(roii)).^2 + ...
      (YGRID - roidata.centery(roii)).^2) ...
      <= roidata.radii(roii)^2 );
  end
  
  if didtrackwings,
    [wingtrx,wingperframedata,wingtrackinfo,wingperframeunits] = TrackWingsHelper(trackdata.trx,moviefile,double(bgmed),isarena,params.wingtracking_params,...
      'firstframe',params.firstframetrack,...
      'debug',params.DEBUG,...
      'framestrack',unique(cat(2,framestrack{:})),...
      'perframedata',trackdata.perframedata);
  else  
    [wingtrx,wingperframedata,wingtrackinfo,wingperframeunits] = TrackWingsHelper(trackdata.trx,moviefile,double(bgmed),isarena,params.wingtracking_params,...
      'firstframe',params.firstframetrack,...
      'debug',params.DEBUG);
  end
  trackdata.trackwings_timestamp2 = wingtrackinfo.trackwings_timestamp;
  trackdata.trackwings_version = wingtrackinfo.trackwings_version;
  trackdata.trx = wingtrx;
  trackdata.perframedata = wingperframedata;
  trackdata.perframeunits = wingperframeunits;

end

end

%% convert to real units

if isfield(roidata,'pxpermm'),
  
  dorotate = isfield(roidata,'rotateby');
  dotranslate = all(isfield(roidata,{'centerx','centery'}));
  if dorotate,
    costheta = cos(roidata.rotateby);
    sintheta = sin(roidata.rotateby);
    R = [costheta,-sintheta;sintheta,costheta];
  end
  
  for fly = 1:numel(trackdata.trx),
    
    roii = trackdata.trx(fly).roi;
    x = trackdata.trx(fly).x;
    if dotranslate,
      x = x - roidata.centerx(roii);
    end
    x = x / roidata.pxpermm;
    y = trackdata.trx(fly).y;
    if dotranslate,
      y = y - roidata.centery(roii);
    end
    y = y / roidata.pxpermm;
    a = trackdata.trx(fly).a / roidata.pxpermm;
    b = trackdata.trx(fly).b / roidata.pxpermm;
    theta = trackdata.trx(fly).theta;
    if dorotate,
      p = R*[x;y];
      x = p(1,:);
      y = p(2,:);
      theta = modrange(theta + roidata.rotateby,-pi,pi);
    end
      
    trackdata.trx(fly).x_mm = x;
    trackdata.trx(fly).y_mm = y;
    trackdata.trx(fly).theta_mm = theta;
    trackdata.trx(fly).a_mm = a;
    trackdata.trx(fly).b_mm = b;
    
    trackdata.trx(fly).pxpermm = roidata.pxpermm;
    trackdata.trx(fly).fps = 1/median(timestamps);
        
  end
  
end

dt = diff(timestamps);
if all(isnan(dt)),
  fps = 30;
  warning('Unknown fps, assigning to 30');
  mediandt = 1/fps;
else
  mediandt = median(dt(~isnan(dt)));
  fps = 1/mediandt;
end
for fly = 1:nflies,
  trackdata.trx(fly).fps = fps;
end

if isfield(params,'usemediandt') && params.usemediandt,
  
  for fly = 1:nflies,
    trackdata.trx(fly).dt = repmat(mediandt,[1,trackdata.trx(fly).nframes-1]);
  end
  
end


%% add arena parameters

if all(isfield(roidata,{'centerx','centery','radii'})),

  for i = 1:nflies,
    roii = trackdata.trx(i).roi;
    trackdata.trx(i).arena.arena_radius_mm = roidata.radii(roii) / roidata.pxpermm;
    trackdata.trx(i).arena.arena_center_mm_x = 0;
    trackdata.trx(i).arena.arena_center_mm_y = 0;
  end
  
end
    
%% clean up

fprintf('Clean up...\n');

if fid > 1,
  try
    fclose(fid);
  catch ME,
    warning('Could not close movie: %s',getReport(ME));
  end
end

if exist(tmpfilename,'file'),
  try
    delete(tmpfilename);
  catch ME,
    warning('Could not delete tmp file: %s',getReport(ME));
  end
end