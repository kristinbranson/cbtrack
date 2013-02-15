function trackdata = CourtshipBowlTrack(expdir,varargin)

version = '0.1';
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
trackdata = rmfield(trackdata,{'trx','timestamps'});
outfilename = fullfile(expdir,dataloc_params.trackingdatamatfilestr);
save(outfilename,'-struct','trackdata');
outfilename = fullfile(expdir,dataloc_params.trxfilestr);
save(outfilename,'trx','timestamps');
