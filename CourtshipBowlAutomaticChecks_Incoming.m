function [success,msgs,iserror] = CourtshipBowlAutomaticChecks_Incoming(expdir,varargin)

version = '0.11';
timestamp = datestr(now,TimestampFormat);

ParseCourtshipBowlParams;

success = true;
msgs = {};

[DEBUG] = ...
  myparse(leftovers,...
  'debug',false);

check_params = cbparams.auto_checks_incoming;

%% open log file

if isfield(cbparams.dataloc,'automaticchecks_incoming_log') && ~DEBUG,
  logfile = fullfile(expdir,cbparams.dataloc.automaticchecks_incoming_log.filestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

fprintf(logfid,'\n\n***\nRunning CourtshipBowlAutomaticChecks_Incoming version %s analysis_protocol %s at %s\n',version,analysis_protocol,timestamp);

%%

try


metadatafile = fullfile(expdir,cbparams.dataloc.metadata.filestr);
moviefile = fullfile(expdir,cbparams.dataloc.movie.filestr);
temperaturefile = fullfile(expdir,cbparams.dataloc.temperature.filestr); %#ok<NASGU>
outfile = fullfile(expdir,cbparams.dataloc.automaticchecksincomingresults.filestr);

% order matters here: higher up categories have higher priority
categories = {'flag_aborted_set_to_1',...
  'missing_video',...
  'missing_metadata_file',...
  'missing_metadata_fields',...
  'short_video',...
  'missing_capture_files',...
  'flag_redo_set_to_1',...
  'fliesloaded_time_too_short',...
  'fliesloaded_time_too_long',...
  'shiftflytemp_time_too_long',...
  'no_barcode',...
  'incoming_checks_other'};
category2idx = struct;
for i = 1:numel(categories),
  category2idx.(categories{i}) = i;
end
iserror = false(1,numel(categories));

%% check for notstarted
[~,expname] = fileparts(expdir);
if ~isempty(regexp(expname,'notstarted','once')),
  iserror(category2idx.missing_video) = true;
  success = false;
  msgs{end+1} = 'Capture not started';
end

%% read metadata

ismetadata = exist(metadatafile,'file');
if ~ismetadata,
  success = false;
  msgs{end+1} = 'Missing Metadata.xml file';
  iserror(category2idx.missing_metadata_file) = true;
  metadata = struct;
else
  try
    metadata = ReadMetadataFile(metadatafile);
  catch %#ok<CTCH>
    msgs{end+1} = 'Error reading Metadata file';
    success = false;
  end
end
  

%% check for metadata fields

required_fns = {'flag_aborted','flag_redo','seconds_fliesloaded',...
  'screen_type','screen_reason'};
if isfield(check_params,'required_fns'),
  required_fns = check_params.required_metadata_fields;
end
ismissingfn = ~ismember(required_fns,fieldnames(metadata));
if any(ismissingfn),
  success = false;
  msgs{end+1} = ['Missing metadata fields:',sprintf(' %s',required_fns{ismissingfn})];
  iserror(category2idx.missing_metadata_fields) = true;
end

%% check for flags

if isfield(metadata,'flag_aborted') && metadata.flag_aborted ~= 0,
  success = false;
  msgs{end+1} = 'Experiment aborted.';
  iserror(category2idx.flag_aborted_set_to_1) = true;
end

if isfield(metadata,'flag_redo') && metadata.flag_redo ~= 0,
  success = false;
  msgs{end+1} = 'Redo flag set to 1.';
  iserror(category2idx.flag_redo_set_to_1) = true;
end

%% check loading time

if isfield(metadata,'seconds_fliesloaded'),
  if metadata.seconds_fliesloaded < check_params.min_seconds_fliesloaded,
    success = false;
    msgs{end+1} = sprintf('Load time = %f < %f seconds.',metadata.seconds_fliesloaded,check_params.min_seconds_fliesloaded);
    iserror(category2idx.fliesloaded_time_too_short) = true;
  end
  if metadata.seconds_fliesloaded > check_params.max_seconds_fliesloaded,
    success = false;
    msgs{end+1} = sprintf('Load time = %f > %f seconds.',metadata.seconds_fliesloaded,check_params.max_seconds_fliesloaded);
    iserror(category2idx.fliesloaded_time_too_long) = true;
  end
end

%% check video length

if ~exist(moviefile,'file'),
  success = false;
  msgs{end+1} = 'Missing movie file.';
  iserror(category2idx.missing_capture_files) = true;
else
  try
    headerinfo = ufmf_read_header(fullfile(expdir,cbparams.dataloc.movie.filestr));
    nframes = headerinfo.nframes;
    fclose(headerinfo.fid);
    if nframes < check_params.min_ufmf_nframes,
      success = false;
      msgs{end+1} = sprintf('Video contains %d < %d frames.',nframes,check_params.min_ufmf_nframes);
      iserror(category2idx.short_video) = true;
    end
  catch ME,
    success = false;
    msgs{end+1} = sprintf('Error reading movie filefile: %s',getReport(ME));
    iserror(category2idx.incoming_checks_other) = true;
  end
end


%% check for missing files

% movie file
fn = cbparams.dataloc.movie.filestr;
isfile = exist(fullfile(expdir,fn),'file');
if ~isfile,
  iserror(category2idx.missing_video) = true;
end

% all other required data capture files
files = fieldnames(cbparams.dataloc);
isrequireddatacap = structfun(@(x) isstruct(x) && isfield(x,'essential') && (x.essential>=1) && isfield(x,'type') && strcmp(x.type,'data_capture'),cbparams.dataloc);
required_files = setdiff(files(isrequireddatacap),'movie');

for i = 1:numel(required_files),
  fn = cbparams.dataloc.(required_files{i}).filestr;
  if any(fn == '*'),
    isfile = ~isempty(dir(fullfile(expdir,fn)));
  else
    isfile = exist(fullfile(expdir,fn),'file');
  end
  if ~isfile,
    msgs{end+1} = sprintf('Missing file %s',fn); %#ok<AGROW>
    success = false;
    iserror(category2idx.missing_capture_files) = true;
  end
end

%% check for missing desired files

isdesireddatacap = structfun(@(x) isstruct(x) && isfield(x,'essential') && (x.essential==0) && isfield(x,'type') && strcmp(x.type,'data_capture'),cbparams.dataloc);
desired_files = files(isdesireddatacap);

for i = 1:numel(desired_files),
  fn = cbparams.dataloc.(desired_files{i}).filestr;
  if any(fn == '*'),
    isfile = ~isempty(dir(fullfile(expdir,fn)));
  else
    isfile = exist(fullfile(expdir,fn),'file');
  end
  if ~isfile,
    msgs{end+1} = sprintf('Missing desired file %s',fn); %#ok<AGROW>
  end
end

%% output results to file

if DEBUG,
  fid = 1;
else
  if exist(outfile,'file'),
    try
      delete(outfile);
    catch ME,
      warning('Could not delete file %s:\n %s',outfile,getReport(ME));
    end
  end
  fid = fopen(outfile,'w');
end
if fid < 0,
  warning('Could not open automatic checks results file %s for writing, just printing to stdout.',outfile);
  fid = 1;
end
if success,
  fprintf(fid,'automated_pf,U\n');
else
  fprintf(fid,'automated_pf,F\n');
  fprintf(fid,'notes_curation,');
  s = sprintf('%s\\n',msgs{:});
  s = s(1:end-2);
  fprintf(fid,'%s\n',s);
  i = find(iserror,1);
  if isempty(i),
    s = 'incoming_checks_other';
  else
    s = categories{i};
  end
  fprintf(fid,'automated_pf_category,%s\n',s);
  
end
% version info
fprintf(fid,'cbautochecksincoming_version,%s\n',version);
fprintf(fid,'cbautochecksincoming_timestamp,%s\n',timestamp);
fprintf(fid,'analysis_protocol,%s\n',real_analysis_protocol);


if ~DEBUG && fid > 1,
  fclose(fid);
end

catch ME,
  success = false;
  msgs = {getReport(ME)};
end
  
%% print results to log file

fprintf(logfid,'Finished running CourtshipBowlAutomaticChecks_Incoming at %s.\n',datestr(now,'yyyymmddTHHMMSS'));
fprintf(logfid,'success = %d\n',success);
if isempty(msgs),
  fprintf(logfid,'No error or warning messages.\n');
else
  fprintf(logfid,'Warning/error messages:\n');
  fprintf(logfid,'%s\n',msgs{:});
end

if logfid > 1,
  fclose(logfid);
end