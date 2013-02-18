function trackdata = CourtshipBowlTrack(expdir,varargin)

version = '0.1.1';
timestamp = datestr(now,TimestampFormat);

%% parse inputs
ParseCourtshipBowlParams;

% read parameters
params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.trackingparamsfilestr));
SetBackgroundTypes;
if ischar(params.bgmode) && isfield(bgtypes,params.bgmode),
  params.bgmode = bgtypes.(params.bgmode);
end
restart = '';
for i = 1:2:numel(leftovers)-1, %#ok<USENS>
  if strcmpi(leftovers{i},'restart'),
    restart = leftovers{i+1};
  else
    if isdeployed && ~isempty(leftovers{i+1}) && ischar(leftovers{i+1}),
      tmp = str2double(leftovers{i+1});
    else
      tmp = leftovers{i+1};
    end
    params.(leftovers{i}) = tmp;
  end
end
if ~isfield(params,'DEBUG'),
  params.DEBUG = 0;
end
if params.dotrackwings || strcmp(params.assignidsby,'wingsize'),
  params.wingtracking_params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.wingtrackingparamsfilestr));
end

%% load background model

try
  load(fullfile(expdir,dataloc_params.bgmatfilestr),'bgmed');
catch ME,
  error('Could not load background model from file: %s',getReport(ME));
end

%% load roi info

try
  roidata = load(fullfile(expdir,dataloc_params.roidatamatfilestr));
  roidata.nrois = numel(roidata.centerx);
catch ME,
  error('Could not load roi data from file: %s',getReport(ME));
end

%% main function

moviefile = fullfile(expdir,dataloc_params.moviefilestr);
tmpfilename = fullfile(expdir,sprintf('TmpResultsTrackTwoFlies_%s.mat',datestr(now,'yyyymmddTHHMMSSPFFF')));
trackdata = TrackTwoFlies(moviefile,bgmed,roidata,params,'restart',restart,'tmpfilename',tmpfilename);

%% save results

trackdata.courtshipbowltrack_version = version;
trackdata.courtshipbowltrack_timestamp = timestamp;
trackdata.params = params;
trx = trackdata.trx; %#ok<NASGU>
timestamps = trackdata.timestamps; %#ok<NASGU>

% trx
outfilename = fullfile(expdir,dataloc_params.trxfilestr);
if exist(outfilename,'file'),
  delete(outfilename);
end
save(outfilename,'trx','timestamps');

% perframe data
perframedir = fullfile(expdir,dataloc_params.perframedir);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end
perframefns = fieldnames(trackdata.perframedata);
for i = 1:numel(perframefns),
  perframefn = perframefns{i};
  filename = fullfile(perframedir,[perframefn,'.mat']);
  if exist(filename,'file'),
    delete(filename);
  end
  data = trackdata.perframedata.(perframefn); %#ok<NASGU>
  units = trackdata.perframeunits.(perframefn); %#ok<NASGU>
  save(filename,'data','units');
end
% also save sex
perframefns = {'sex','x_mm','y_mm','a_mm','b_mm','theta_mm','x','y','a','b','theta','timestamps','dt'};
for i = 1:numel(perframefns),
  perframefn = perframefns{i};
  filename = fullfile(perframedir,[perframefn,'.mat']);
  if strcmp(perframefn,'sex') && ~isfield(trackdata.trx,'sex'),
    data = cell(1,numel(trackdata.trx));
    for fly = 1:numel(trackdata.trx),
      data{fly} = repmat('?',[1,trackdata.trx(fly).nframes]);
    end
    units = parseunits('unit'); %#ok<NASGU>
  elseif ~isfield(trackdata.trx,perframefn) || ~isfield(trackdata.perframeunits,perframefn),
    continue;
  else
    data = {trackdata.trx.(perframefn)}; %#ok<NASGU>
    units = trackdata.perframeunits.(perframefn);     %#ok<NASGU>
  end
  save(filename,'data','units');
  
end


% tracking data without the trx
trackdata = rmfield(trackdata,{'trx','timestamps','perframedata','perframeunits'});
outfilename = fullfile(expdir,dataloc_params.trackingdatamatfilestr);
save(outfilename,'-struct','trackdata');


