% ScriptHistogramRelativePositions20130502

%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;

expfile = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/expdirlists/grouped_expdirs_galit_rejection_wingclipped_20130207.txt';
analysis_protocol = '20130212_galit_rejection_wingclipped_20130207';

%% read in the experiment list

expdirs = ReadGroupedExperimentList(expfile);

%% histogram

conditions = fieldnames(expdirs);
allmeannormfrac = struct;
for i = 1:numel(conditions),
  condition = conditions{i};
  expdirscurr = expdirs.(condition);
  for j = 1:numel(expdirscurr),
    [meannormfrac_curr,~,~,~,edges_angle,edges_dist] = ...
      HistogramRelativePositions(expdirscurr{j});
    if j == 1,
      meannormfrac = meannormfrac_curr;
      fns = fieldnames(meannormfrac);
    else
      for k = 1:numel(fns),
        meannormfrac.(fns{k}) = meannormfrac.(fns{k}) + meannormfrac_curr.(fns{k});
      end
    end
  end
  for k = 1:numel(fns),
    meannormfrac.(fns{k}) = meannormfrac.(fns{k}) / numel(expdirscurr);
  end
  allmeannormfrac.(condition) = meannormfrac;
end

%% plot

figure(3);
clf;
hax = createsubplots(2,numel(conditions),[.02,.025;0,.01]);
hax = reshape(hax,[2,numel(conditions)]);
for i = 1:numel(conditions),
  condition = conditions{i};
  fns = fieldnames(allmeannormfrac.(condition));
  for j = 1:numel(fns),
    axes(hax(j,i));
    polarimagesc(edges_dist,edges_angle+pi/2,[0,0],allmeannormfrac.(condition).(fns{j}));
    axis equal;
    axisalmosttight;
    title(sprintf('%s, %s',condition,fns{j}),'Interpreter','none');
    colorbar;
  end
end
set(hax(1,:),'XTickLabel',{});
set(hax(:,2:end),'YTickLabel',{});

%% histogram distance only

conditions = fieldnames(expdirs);
allmeanfrac_dist = struct;
allstdfrac_dist = struct;
allstderrfrac_dist = struct;
allfracperfly_dist = struct;
allfly2exp = struct;
for i = 1:numel(conditions),
  condition = conditions{i};
  expdirscurr = expdirs.(condition);
  meanfrac_dist = 0;
  fracperfly = [];
  fly2exp = [];
  meanfrac_dist_perexp = struct;
  for j = 1:numel(expdirscurr),
    [meanfrac_curr,fracperfly_curr,centers_dist] = ...
      HistogramDistance(expdirscurr{j});
    meanfrac_dist = meanfrac_dist + meanfrac_curr;
    nflies_curr = size(fracperfly_curr,2);
    fracperfly(:,end+1:end+nflies_curr) = fracperfly_curr;
    fly2exp(end+1:end+nflies_curr) = j;
  end
  meanfrac_dist = meanfrac_dist / numel(expdirscurr);
  stdfrac_dist = std(fracperfly,1,2);
  stderrfrac_dist = stdfrac_dist / sqrt(size(fracperfly,2));
  allmeanfrac_dist.(condition) = meanfrac_dist;
  allstdfrac_dist.(condition) = stdfrac_dist;
  allstderrfrac_dist.(condition) = stderrfrac_dist;
  allfracperfly_dist.(condition) = fracperfly;
  allfly2exp.(condition) = fly2exp;
end

%% plot

figure(4);
clf;
hold on;
colors = lines(numel(conditions));

% % plot individual flies
% for i = 1:numel(conditions),
%   condition = conditions{i};
%   nexpscurr = numel(expdirs.(condition));
%   tmp = linspace(.25,.5,nexpscurr)';
%   colorscurr = bsxfun(@plus,bsxfun(@times,colors(i,:),1-tmp),tmp);
%   for j = 1:nexpscurr,
%     plot(centers_dist,allfracperfly_dist.(condition)(:,allfly2exp.(condition)==j),...
%       '.-','Color',colorscurr(j,:));
%   end
% end

% plot standard errors
for i = 1:numel(conditions),
  condition = conditions{i};
  errorpatch(centers_dist,allmeanfrac_dist.(condition),...
    allstderrfrac_dist.(condition),'FaceColor',...
    colors(i,:)*.5+.5);
end
% plot condition means
h = nan(1,numel(conditions));
for i = 1:numel(conditions),
  condition = conditions{i};
  h(i) = plot(centers_dist,allmeanfrac_dist.(condition),'.-',...
    'Color',colors(i,:),'LineWidth',2);
end  
axisalmosttight;
legend(h,conditions,'Interpreter','none','Location','NorthWest');
xlabel('Distance between fly centroids (mm)');
ylabel('Mean fraction of time');
title('Mean and standard error of histogram of inter-fly distance');
