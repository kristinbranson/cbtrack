%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;

expfile = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/expdirlists/grouped_expdirs_galit_rejection_wingclipped_20130207.txt';
analysis_protocol = '20130212_galit_rejection_wingclipped_20130207';

%% read in the experiment list

expdirs = ReadGroupedExperimentList(expfile);

%% auto checks

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    if i == 1 && j == 2, continue; end
    expdir = expdirs.(fn){j};
    CourtshipBowlAutomaticChecks_Incoming(expdir,'analysis_protocol',analysis_protocol);
  end
end

%% detect arenas

fns = fieldnames(expdirs);
roidata = struct;
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    roidata.(fn){j} = CourtshipBowlDetectROIs(expdir,'analysis_protocol',analysis_protocol);
  end
end

%% track

trackdata = CourtshipBowlTrack(expdir,'analysis_protocol',analysis_protocol,'DEBUG',0);

%% per-frame features

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    if i == 1 && j == 2, continue; end
    expdir = expdirs.(fn){j};
    CourtshipBowlComputePerFrameFeatures(expdir,'analysis_protocol',analysis_protocol,'forcecompute',false);
  end
end

%% results movie

fns = fieldnames(expdirs);
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    if i == 1 && j == 2, continue; end
    expdir = expdirs.(fn){j};
    disp(expdir);
    CourtshipBowlMakeResultsMovie(expdir,'analysis_protocol',analysis_protocol);
  end
end
