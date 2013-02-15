%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;

expfile = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/expdirs_galit_rejection_wingclipped_20130207.txt';
analysis_protocol = '20130212_galit_rejection_wingclipped_20130207';

%% read in the experiment list

expdirs = ReadGroupedExperimentList(expfile);

%% detect arenas

roidata = struct;
for i = 1:numel(fns),
  fn = fns{i};
  for j = 1:numel(expdirs.(fn)),
    expdir = expdirs.(fn){j};
    roidata.(fn){j} = CourtshipBowlDetectROIs(expdir,'analysis_protocol',analysis_protocol);
  end
end

%% track

trackdata = CourtshipBowlTrack(expdir,'analysis_protocol',analysis_protocol,'DEBUG',1);
