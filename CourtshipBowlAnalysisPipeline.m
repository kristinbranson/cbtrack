function [success,msgs,stage] = CourtshipBowlAnalysisPipeline(expdir,varargin)

success = false;
msgs = {};
stage = 'start';

ParseCourtshipBowlParams;

[forcecompute,do_auto_checks_incoming,do_detect_rois,do_track,...
  do_compute_perframe_features,do_results_movie,...
  do_auto_checks_complete,...
  auto_checks_incoming_params,...
  detect_rois_params,...
  track_params,...
  compute_perframe_features_params,...
  results_movie_params,...
  auto_checks_complete_params,...
  DEBUG] = ...
  myparse(varargin,...
  'forcecompute',false,...
  'do_auto_checks_incoming',true,...
  'do_detect_rois',true,...
  'do_track',true,...
  'do_compute_perframe_features',true,...
  'do_results_movie',true,...
  'do_auto_checks_complete',true,...
  'auto_checks_incoming_params',{},...
  'detect_rois_params',{},...
  'track_params',{},...
  'compute_perframe_features_params',{},...
  'results_movie_params',{},...
  'auto_checks_complete_params',{},...
  'debug',0);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end
if ischar(do_auto_checks_incoming),
  do_auto_checks_incoming = str2double(do_auto_checks_incoming);
end
if ischar(do_detect_rois),
  do_detect_rois = str2double(do_detect_rois);
end
if ischar(do_track),
  do_track = str2double(do_track);
end
if ischar(do_compute_perframe_features),
  do_compute_perframe_features = str2double(do_compute_perframe_features);
end
if ischar(do_results_movie),
  do_results_movie = str2double(do_results_movie);
end
if ischar(do_auto_checks_complete),
  do_auto_checks_complete = str2double(do_auto_checks_complete);
end

%% check that experiment exists

if ~exist(expdir,'dir'),
  msgs = {sprintf('Experiment directory %s does not exist',expdir)};
  fprintf('%s\n',msgs{:});
  return;
end

%% incoming checks

stage = 'auto_checks_incoming';

if do_auto_checks_incoming,
  [todo_required,~,todo_desired,~] = CheckForMissingFiles(expdir,cbparams,stage);
  todo = todo_required || todo_desired;
  if forcecompute || todo,
    try
      fprintf('AutomaticChecks_Incoming...\n');
      [success1,msgs] = CourtshipBowlAutomaticChecks_Incoming(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'debug',DEBUG,...
        auto_checks_incoming_params{:});
      if ~success1,
        fprintf('AutomaticChecks_Incoming failed:\n');
        fprintf('%s\n',msgs{:});
        return;
      end      
    catch ME,
      msgs = {sprintf('Error running AutomaticChecks_Incoming:\n%s',getReport(ME))};
      fprintf('AutomaticChecks_Incoming failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingrequiredfile,missingrequiredfiles,...
    ismissingdesiredfile,missingdesiredfiles] = CheckForMissingFiles(expdir,cbparams,stage);
  if ismissingrequiredfile,
    msgs = cellfun(@(x) sprintf('Missing required auto_checks_incoming file %s',x),missingrequiredfiles,'UniformOutput',false);
    fprintf('AutomaticChecks_Incoming failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
  if ismissingdesiredfile,
    warning(['Missing desired auto_checks_incoming files: ',sprintf('%s ',missingdesiredfiles{:})]);
  end
  
end

%% detect ROIs

stage = 'detect_rois';

if do_detect_rois,
  [todo_required,~,todo_desired,~] = CheckForMissingFiles(expdir,cbparams,stage);
  todo = todo_required || todo_desired;
  if forcecompute || todo,
    try
      fprintf('DetectROIs...\n');
      roidata = CourtshipBowlDetectROIs(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'debug',DEBUG,...
        detect_rois_params{:});
      drawnow;
      
      fprintf('Found %d rois.\n\n',numel(roidata.nflies_per_roi));
      fprintf('nflies\tnrois\n');
      for i = 0:2,
        fprintf('%d\t%d\n',i,nnz(roidata.nflies_per_roi==i));
      end
      fprintf('>2\t%d\n',nnz(roidata.nflies_per_roi>2));
      fprintf('\n');
      
    catch ME,
      msgs = {sprintf('Error running DetectROIs:\n%s',getReport(ME))};
      fprintf('DetectROIs failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure detect_rois files exist
  [ismissingrequiredfile,missingrequiredfiles,...
    ismissingdesiredfile,missingdesiredfiles] = CheckForMissingFiles(expdir,cbparams,stage);
  if ismissingrequiredfile,
    msgs = cellfun(@(x) sprintf('Missing required detect_rois file %s',x),missingrequiredfiles,'UniformOutput',false);
    fprintf('ComputePerFrameFeatures failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
  if ismissingdesiredfile,
    warning(['Missing desired detect_rois files: ',sprintf('%s ',missingdesiredfiles{:})]);
  end
  
end

%% track

stage = 'track';

if do_track,
  [todo_required,~,todo_desired,~] = CheckForMissingFiles(expdir,cbparams,stage);
  todo = todo_required || todo_desired;
  if forcecompute || todo,
    try
      fprintf('Track...\n');
      CourtshipBowlTrack(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'debug',DEBUG,...
        track_params{:});
      drawnow;
      
    catch ME,
      msgs = {sprintf('Error running Track:\n%s',getReport(ME))};
      fprintf('Track failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure track files exist
  [ismissingrequiredfile,missingrequiredfiles,...
    ismissingdesiredfile,missingdesiredfiles] = CheckForMissingFiles(expdir,cbparams,stage);
  if ismissingrequiredfile,
    msgs = cellfun(@(x) sprintf('Missing required track file %s',x),missingrequiredfiles,'UniformOutput',false);
    fprintf('ComputePerFrameFeatures failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
  if ismissingdesiredfile,
    warning(['Missing desired track files: ',sprintf('%s ',missingdesiredfiles{:})]);
  end
  
end

%% per-frame features

stage = 'compute_perframe_features';

if do_compute_perframe_features,
  [todo_required,~,todo_desired,~] = CheckForMissingFiles(expdir,cbparams,stage);
  todo = todo_required || todo_desired;
  if forcecompute || todo,
    try
      fprintf('Track...\n');
      CourtshipBowlComputePerFrameFeatures(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'debug',DEBUG,...
        compute_perframe_features_params{:});
      drawnow;
      
    catch ME,
      msgs = {sprintf('Error running ComputePerFrameFeatures:\n%s',getReport(ME))};
      fprintf('ComputePerFrameFeatures failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure compute_perframe_features files exist
  [ismissingrequiredfile,missingrequiredfiles,...
    ismissingdesiredfile,missingdesiredfiles] = CheckForMissingFiles(expdir,cbparams,stage);
  if ismissingrequiredfile,
    msgs = cellfun(@(x) sprintf('Missing required compute_perframe_features file %s',x),missingrequiredfiles,'UniformOutput',false);
    fprintf('ComputePerFrameFeatures failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
  if ismissingdesiredfile,
    warning(['Missing desired compute_perframe_features files: ',sprintf('%s ',missingdesiredfiles{:})]);
  end
  
end

%% results movie

stage = 'results_movie';

if do_results_movie,
  [todo_required,~,todo_desired,~] = CheckForMissingFiles(expdir,cbparams,stage);
  todo = todo_required || todo_desired;
  if forcecompute || todo,
    try
      fprintf('AutomaticChecks_Incoming...\n');
      CourtshipBowlMakeResultsMovie(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'debug',DEBUG,...
        results_movie_params{:});
    catch ME,
      msgs = {sprintf('Error running MakeResultsMovie:\n%s',getReport(ME))};
      fprintf('MakeResultsMovie failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure results movie files exist
  [ismissingrequiredfile,missingrequiredfiles,...
    ismissingdesiredfile,missingdesiredfiles] = CheckForMissingFiles(expdir,cbparams,stage);
  if ismissingrequiredfile,
    msgs = cellfun(@(x) sprintf('Missing required results_movie file %s',x),missingrequiredfiles,'UniformOutput',false);
    fprintf('MakeResultsMovie failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
  if ismissingdesiredfile,
    warning(['Missing desired results_movie files: ',sprintf('%s ',missingdesiredfiles{:})]);
  end
  
end

%% incoming checks

stage = 'auto_checks_complete';

if do_auto_checks_complete,
  [todo_required,~,todo_desired,~] = CheckForMissingFiles(expdir,cbparams,stage);
  todo = todo_required || todo_desired;
  if forcecompute || todo,
    try
      fprintf('AutomaticChecks_Complete...\n');
      [success1,msgs] = CourtshipBowlAutomaticChecks_Complete(expdir,...
        'settingsdir',settingsdir,'analysis_protocol',analysis_protocol,...
        'debug',DEBUG,...
        auto_checks_complete_params{:});
      if ~success1,
        fprintf('AutomaticChecks_Complete failed:\n');
        fprintf('%s\n',msgs{:});
        return;
      end      
    catch ME,
      msgs = {sprintf('Error running AutomaticChecks_Complete:\n%s',getReport(ME))};
      fprintf('AutomaticChecks_Complete failed:\n');
      fprintf('%s\n',msgs{:});
      return;
    end
  end
  
  % make sure automatic checks files exist
  [ismissingrequiredfile,missingrequiredfiles,...
    ismissingdesiredfile,missingdesiredfiles] = CheckForMissingFiles(expdir,cbparams,stage);
  if ismissingrequiredfile,
    msgs = cellfun(@(x) sprintf('Missing required auto_checks_complete file %s',x),missingrequiredfiles,'UniformOutput',false);
    fprintf('AutomaticChecks_Complete failed:\n');
    fprintf('%s\n',msgs{:});
    return;
  end  
  if ismissingdesiredfile,
    warning(['Missing desired auto_checks_complete files: ',sprintf('%s ',missingdesiredfiles{:})]);
  end
  
end

%% done

success = true;
