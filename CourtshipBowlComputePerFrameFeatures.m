function pfdata = CourtshipBowlComputePerFrameFeatures(expdir,varargin)

version = '0.1';
timestamp = datestr(now,TimestampFormat);

pfdata = struct;
pfdata.cbdetectrois_version = version;
pfdata.cbdetectrois_timestamp = timestamp;

%% parameters

ParseCourtshipBowlParams;
[forcecompute,perframefns,DEBUG] = ...
  myparse(leftovers,'forcecompute',true,'perframefns',{},'DEBUG',false);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end
if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
end

perframe_params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.perframeparamsfilestr));

%% load the trx

fprintf('Initializing trx...\n');

% uses the Trx code in JAABA
trx = Trx('trxfilestr',dataloc_params.trxfilestr,...
  'perframedir',dataloc_params.perframedir,...
  'moviefilestr',dataloc_params.moviefilestr,...
  'perframe_params',perframe_params);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir,'openmovie',false);

%% log file

if isfield(dataloc_params,'perframefeature_logfilestr') && ~DEBUG,
  logfile = fullfile(expdir,dataloc_params.perframefeature_logfilestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf(logfid,'\n\n***\nRunning FlyBowlComputePerFrameFeatures version %s analysis_protocol %s at %s\n',version,analysis_protocol,timestamp);


%% compute per-frame features

if isempty(perframefns) 
  perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
  perframefns = importdata(perframefnsfile);
end
nfns = numel(perframefns);

% clean this data to force computation
if forcecompute,
  %deletefns = setdiff(perframefns,Trx.TrajectoryFieldNames());
  trx.CleanPerFrameData();
end

% compute each
for i = 1:nfns,
  fn = perframefns{i};
  fprintf(logfid,'Computing %s...\n',fn);
  trx.(fn); 
end

%% close log

fprintf(logfid,'Finished running FlyBowlComputePerFrameFeatures at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end