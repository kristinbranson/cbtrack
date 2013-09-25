function trackdata = CourtshipBowlTrack(expdir,varargin)

version = '0.1.2';
timestamp = datestr(now,TimestampFormat);

%% parse inputs
ParseCourtshipBowlParams;
params = cbparams.track;
metadatafile = fullfile(expdir,cbparams.dataloc.metadata.filestr);
metadata = ReadMetadataFile(metadatafile);

SetBackgroundTypes;
if ischar(params.bgmode) && isfield(bgtypes,params.bgmode),
  params.bgmode = bgtypes.(params.bgmode);
end
restart = '';
for i = 1:2:numel(leftovers)-1, %#ok<NODEF>
  if strcmpi(leftovers{i},'restart'),
    restart = leftovers{i+1};
  else
    
    if strcmpi(leftovers{i},'debug'),
      leftovers{i} = 'DEBUG'; %#ok<AGROW>
    end
    
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
  params.wingtracking_params = cbparams.wingtrack;
end

%% open log file

if isfield(cbparams.dataloc,'track_log'),
  logfile = fullfile(expdir,cbparams.dataloc.track_log.filestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

fprintf(logfid,'\n\n***\nRunning CourtshipBowlTrack version %s analysis_protocol %s at %s\n',version,real_analysis_protocol,timestamp);

%% load background model

fprintf(logfid,'Loading background model...\n');

try
  load(fullfile(expdir,cbparams.dataloc.bgmat.filestr),'bgmed');
catch ME,
  fprintf(logfid,'Could not load background model from file: %s',getReport(ME));
  error('Could not load background model from file: %s',getReport(ME));
end

%% load roi info

fprintf(logfid,'Loading ROI info...\n');

try
  roidata = load(fullfile(expdir,cbparams.dataloc.roidatamat.filestr));
  roidata.nrois = numel(roidata.centerx);
catch ME,
  fprintf(logfid,'Could not load roi data from file: %s',getReport(ME));
  error('Could not load roi data from file: %s',getReport(ME));
end

%% main function

fprintf(logfid,'Calling TrackTwoFlies...\n');
moviefile = fullfile(expdir,cbparams.dataloc.movie.filestr);
tmpfilename = fullfile(expdir,sprintf('TmpResultsTrackTwoFlies_%s.mat',datestr(now,'yyyymmddTHHMMSSPFFF')));
trackdata = TrackTwoFlies(moviefile,bgmed,roidata,params,'restart',restart,'tmpfilename',tmpfilename,'logfid',logfid);

%% save results

trackdata.courtshipbowltrack_version = version;
trackdata.courtshipbowltrack_timestamp = timestamp;
trackdata.analysis_protocol = real_analysis_protocol;
trackdata.params = params;
trx = trackdata.trx; %#ok<NASGU>
timestamps = trackdata.timestamps; %#ok<NASGU>

% trx
outfilename = fullfile(expdir,cbparams.dataloc.trx.filestr);
fprintf(logfid,'Saving trx to file %s...\n',outfilename);
if exist(outfilename,'file'),
  delete(outfilename);
end
save(outfilename,'trx','timestamps');

% perframe data
perframedir = fullfile(expdir,cbparams.dataloc.perframedir.filestr);
fprintf(logfid,'Saving a bit of per-frame data to directory %s...\n',perframedir);
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
    if isfield(metadata,'gender') && ismember(lower(metadata.gender),{'m','b'}),
      gender = metadata.gender;
    else
      gender = '?';
    end
    for fly = 1:numel(trackdata.trx),
      data{fly} = repmat(gender,[1,trackdata.trx(fly).nframes]);
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
outfilename = fullfile(expdir,cbparams.dataloc.trackingdatamat.filestr);
save(outfilename,'-struct','trackdata');

%% close log file

fprintf(logfid,'Finished running CourtshipBowlTrack at %s.\n',datestr(now,'yyyymmddTHHMMSS'));
if logfid > 1,
  fclose(logfid);
end
