if ~exist(expdir,'dir'),
  error('Experiment directory %s does not exist',expdir);
end

[analysis_protocol,settingsdir,paramsfilestr,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/settings',...
  'paramsfilestr','params.xml');

analysis_protocol_dir = fullfile(settingsdir,analysis_protocol);
real_analysis_protocol_dir = readunixlinks(analysis_protocol_dir);
[~,real_analysis_protocol] = fileparts(analysis_protocol_dir);
cbparams = ReadXMLParams(fullfile(analysis_protocol_dir,paramsfilestr));
