%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;

expfile = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/expdirlists/expdirs_shelby_cbtest_20130715.txt';
analysis_protocol = '20130731_shelby';
analysis_protocol_base = '20130604';
rootdatadirs = {'/groups/branson/bransonlab/projects/CourtshipBowls/data/Shelby/10minVideos'
  '/groups/branson/bransonlab/projects/CourtshipBowls/data/Shelby/1hrVideos'};

logfid = 1;

%% data

expdirs_pertype = ReadGroupedExperimentList(expfile);
expdirs = {};
fns = fieldnames(expdirs_pertype);
for i = 1:numel(fns),
  expdirs = [expdirs,expdirs_pertype.(fns{i})];
end

%% parameters

% for finding rois
maxroidcenter = 25;
roiradiusfactor = .975;

%% roi detection parameters

roidata_all = {};

varargin = {'analysis_protocol',analysis_protocol_base};
ParseCourtshipBowlParams;
roiparams = cbparams.detect_rois;
tracking_params = cbparams.track;
roiparams.roimus = [roiparams.roimus.x(:),roiparams.roimus.y(:)];

%% loop over experiments

for expi = 1:numel(expdirs),
  
  expdir = expdirs{expi};
  [~,experiment_name] = fileparts(expdir);
  
  %% try using previous detection parameters
  
  % open movie
  moviefile = fullfile(expdir,cbparams.dataloc.movie.filestr);
  fprintf(logfid,'Opening movie file %s...\n',moviefile);
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
  im = readframe(1);
  [imheight,imwidth,~] = size(im);
  
  % estimate the background model
  
  % compute background model
  fprintf(logfid,'Computing background model for %s...\n',experiment_name);
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
  fclose(fid);
  
  % detect rois
  fprintf(logfid,'Detecting rois for %s...\n',experiment_name);
  [roidata.centerx,roidata.centery,roidata.radii,roidata.scores] = DetectCourtshipBowlROIs(bgmed,roiparams);
  nrois = numel(roidata.centerx);
  fprintf(logfid,'Detected %d ROIs with mean radius %f (std = %f, min = %f, max = %f)\n',...
    nrois,mean(roidata.radii),std(roidata.radii,1),min(roidata.radii),max(roidata.radii));
  
  % plot results
  hfig = 18159;
  figure(hfig);
  clf;
  figpos = [10,10,imwidth,imheight];
  set(hfig,'Units','pixel','Position',figpos);
  axes('Position',[0,0,1,1]);
  imagesc(bgmed,[0,255]);
  axis image;
  axis off;
  hold on;
  %colors = jet(nrois)*.7;
  for i = 1:nrois,
    drawellipse(roidata.centerx(i),roidata.centery(i),0,roidata.radii(i),roidata.radii(i),'Color','w');
  end
  
  set(gca,'CLim',[prctile(bgmed(:),2),max(bgmed(:))]);
  colormap jet;
  set(hfig,'Units','pixels','Position',figpos);
  
  %% ask if these are ok
  
  res = questdlg('Use these detections or draw circles manually?','Automatic detections ok?','Use these','Draw manually','Cancel','Use these');
  
  if strcmpi(res,'Cancel'),
    break;
  elseif strcmpi(res,'Use these'),
    roidata_all{expi} = roidata;
  elseif strcmpi(res,'Draw manually'),
    
    figure(hfig);
    clf;
    hax = gca;
    imagesc(bgmed,[0,255]);
    axis image;
    colormap jet;
    hold on;
    
    title('Click to add points on next circle. Double-click on the polygon to continue.');
    
    roidata = struct('centerx',[],'centery',[],'radii',[]);
    
    while true,
      
      h = impoly(hax,'Closed',true);
      position = wait(h);
      % fit a circle to the clicked points
      [xc,yc,radius] = fit_circle_to_points(position(:,1),position(:,2));
      
      hold on;
      hcircle = nan(1,2);
      hcircle(1) = plot(xc,yc,'bx','LineWidth',2);
      % plot the entire circle
      theta = 0:0.01:2*pi;
      Xfit = radius*cos(theta) + xc;
      Yfit = radius*sin(theta) + yc;
      hcircle(2) = plot(Xfit, Yfit,'b-','LineWidth',2);
      hplot = plot(position(:,1),position(:,2),'r.-');
      
      res = questdlg('Choose next action:','ROI next action','Add another ROI','Finish','Undo','Add another ROI');
      switch res,
        case {'Finish','Add another ROI'},
          roidata.centerx(end+1) = xc;
          roidata.centery(end+1) = yc;
          roidata.radii(end+1) = radius;
          if strcmp(res,'Finish'),
            break;
          end
        case 'Undo',
          if ishandle(hplot),
            delete(hplot);
          end
          delete(hcircle(ishandle(hcircle)));
      end
    end
    
    nrois = numel(roidata.centerx);
    roidata_all{expi} = roidata;
    
  end
  
  % update roi parameter estimates
  
  roimus = zeros(0,2);
  roiradii = [];
  roiidx = false(expi,0);
  roins = [];
  for i = 1:expi,
    for j = 1:numel(roidata_all{i}.centerx),
      idxcurr = find(~roiidx(i,:));
      if ~isempty(idxcurr),
        [d,kk] = min(sqrt( (roimus(idxcurr,1)-roidata_all{i}.centerx(j)).^2 + (roimus(idxcurr,2)-roidata_all{i}.centery(j)).^2 ));
        k = idxcurr(kk);
      else
        d = inf;
      end
      if d > maxroidcenter,
        k = size(roimus,1)+1;
        roins(k) = 0; %#ok<SAGROW>
        roimus(k,:) = 0;
        roiradii(k) = 0; %#ok<SAGROW>
      end
      
      roiidx(i,k) = true;
      roimus(k,:) = (roimus(k,:)*roins(k) + [roidata_all{i}.centerx(j),roidata_all{i}.centery(j)]) / (roins(k)+1);
      roiradii(k) = (roiradii(k)*roins(k) + roidata_all{i}.radii(j)) / (roins(k)+1); %#ok<SAGROW>
      roins(k) = roins(k) + 1; %#ok<SAGROW>
    end
  end
  
  meanroiradius = nanmean(roiradii);
  roiparams.meanroiradius = meanroiradius;
  roiparams.roimus = roimus;
end

fprintf('meanroiradius="%f"\n',roiparams.meanroiradius);

fprintf('<roimus\n');
fprintf('\tx="%f',roiparams.roimus(1,1));
fprintf(',%f',roiparams.roimus(2:end,1));
fprintf('"\n');
fprintf('\ty="%f',roiparams.roimus(1,2));
fprintf(',%f',roiparams.roimus(2:end,2));
fprintf('"\n');
fprintf('\\>\n');

%% set background subtraction thresholds

hfig = 1;
figure(hfig);
clf;
hax = gca;

expi = 0;
fid = 0;

him = nan;
h1 = nan;
h2 = nan;
hti = nan;

bgthresh = tracking_params.bgthresh;
bgthresh_low = tracking_params.bgthresh_low;
%%

done = false;
while true,

  if expi == 0,
    
    while true,
      fprintf('Select movie index:\n');
      for i = 1:numel(expdirs),
        fprintf('%d: %s\n',i,expdirs{i});
      end
      expi = input('');
      if expi < 1 || expi > numel(expdirs) || round(expi) ~= expi,
        continue;
      end
      break;
    end
           
    expdir = expdirs{expi};
    [~,experiment_name] = fileparts(expdir);    
    % open movie
    moviefile = fullfile(expdir,cbparams.dataloc.movie.filestr);
    fprintf(logfid,'Opening movie file %s...\n',moviefile);
    [readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
    imheight = headerinfo.nr;
    imwidth = headerinfo.nc;
    
    % estimate the background model

    % compute background model
    fprintf(logfid,'Computing background model for %s...\n',experiment_name);
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

    f = 1;
    
  end
  
  im = readframe(f);
  
  % subtract off background
  switch tracking_params.bgmode,
    case 'DARKBKGD',
      dbkgd = imsubtract(im,bgmed);
    case 'LIGHTBKGD',
      dbkgd = imsubtract(bgmed,im);
    case 'OTHERBKGD',
      dbkgd = imabsdiff(im,bgmed);
  end

  % threshold
  isfore = dbkgd >= bgthresh;
  b = bwboundaries(isfore);
  rc = zeros(0,2);
  for i = 1:numel(b),
    rc = [rc;nan(1,2);b{i}];
  end

  isforelow = dbkgd >= bgthresh_low;
  b = bwboundaries(isforelow);
  rclow = zeros(0,2);
  for i = 1:numel(b),
    rclow = [rclow;nan(1,2);b{i}];
  end

  
  if ishandle(him),
    set(him,'CData',im);
  else
    hold(hax,'off');
    him = imagesc(im,'Parent',hax,[0,255]);
    axis(hax,'image');
    colorbar('Peer',hax);
  end
  if ishandle(h1),
    set(h1,'XData',rclow(:,2),'YData',rclow(:,1));
  else
    hold(hax,'on');
    h1 = plot(hax,rclow(:,2),rclow(:,1),'k.');
  end
  if ishandle(h2),
    set(h2,'XData',rc(:,2),'YData',rc(:,1));
  else
    hold(hax,'on');
    h2 = plot(hax,rc(:,2),rc(:,1),'w.');
  end
  if ishandle(hti),
    set(hti,'String',num2str(f));
  else
    hti = title(hax,num2str(f));
  end
  
  
  while true,
    
    fprintf('1: decrease low threshold\n');
    fprintf('2: increase low threshold\n');
    fprintf('3: decrease high threshold\n');
    fprintf('4: increase high threshold\n');
    fprintf('5: previous frame\n');
    fprintf('6: next frame\n');
    fprintf('7: select frame\n');
    fprintf('8: select random frame\n');
    fprintf('9: select different movie\n');
    fprintf('0: exit\n');
    
    v = input('');
    
    switch v,
      case 1,
        bgthresh_low = bgthresh_low*.99;
      case 2,
        bgthresh_low = bgthresh_high/.99;
      case 3,
        bgthresh = bgthresh*.99;
      case 4,
        bgthresh = bgthresh/.99;
      case 5,
        f = max(1,f-1);
      case 6,
        f = min(nframes,f+1);
      case 7,
        while true,
          f1 = input(sprintf('Enter frame (1 to %d)',nframes));
          if f1 < 1 || f1 > nframes || round(f1) ~= f1,
            continue;
          end
          break;
        end
        f = f1;
      case 9,
        expi = 0;
        if fid > 0,
          fclose(fid);
        end
        fid = 0;
      case 8,
        f = randsample(nframes,1);        
      case 0,
        done = true;
        break;
      otherwise,
        continue;
    end
    
    break;
    
  end

  if done,
    break;
  end
  
end