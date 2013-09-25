%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/simplewing;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/compute_perframe_features;

expfile = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/expdirlists/expdirs_shelby_cbtest_20130715.txt';
%expfile = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/expdirlists/expdirs_cbtest_20121121.txt';
analysis_protocol = '20130731_shelby';

%% read in the experiment list

expdirs = ReadGroupedExperimentList(expfile);

%% auto checks

fns = fieldnames(expdirs);
issuccess = struct;
for i = 1:numel(fns),
  fn = fns{i};
  issuccess.(fn) = false(1,numel(expdirs.(fn)));
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    [success,msgs,iserror] = CourtshipBowlAutomaticChecks_Incoming(expdir,'analysis_protocol',analysis_protocol);
    issuccess.(fn)(j) = success;
    if ~success,
      fprintf('\nAuto checks incoming failed for %s\n',expdir);
      fprintf('%s\n',msgs{:});
    end
      
  end
end

%% detect arenas

fns = fieldnames(expdirs);
roidata = struct;
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    if ~issuccess.(fn)(j),
      continue;
    end
    expdir = expdirs.(fn){j};
    roidata.(fn){j} = CourtshipBowlDetectROIs(expdir,'analysis_protocol',analysis_protocol,'dofixbg',true);
  end
end

%% track

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    CourtshipBowlTrack(expdir,'analysis_protocol',analysis_protocol,'DEBUG',0);
  end
end

%% per-frame features

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    CourtshipBowlComputePerFrameFeatures(expdir,'analysis_protocol',analysis_protocol,'forcecompute',true);
  end
end

%% results movie

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    disp(expdir);
    CourtshipBowlMakeResultsMovie(expdir,'analysis_protocol',analysis_protocol);
  end
end

%% checks

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    disp(expdir);
    [success,msgs,iserror] = CourtshipBowlAutomaticChecks_Complete(expdir,'analysis_protocol',analysis_protocol);
  end
end