function [meanfrac,fracperfly,centers_dist] = HistogramDistance(expdir,varargin)

nbins_dist = 10;
mindist = 1;
maxdist = 25;

[nbins_dist,mindist,maxdist] = myparse(varargin,...
  'nbins_dist',nbins_dist,...
  'mindist',mindist,...
  'maxdist',maxdist);

load(fullfile(expdir,'registered_trx.mat'),'trx');
nflies = numel(trx);
wingtype = cell(1,nflies);
for i = 1:nflies,
  wingtype{i} = trx(i).wingtype{1};
end

[unique_wingtypes,~,wingtypeidx] = unique(wingtype); %#ok<ASGLU>

load(fullfile(expdir,'perframe','dcenter.mat'),'data','units');
dcenter = data;

%edges_dist = logspace(log10(mindist),log10(maxdist),nbins_dist+1);
edges_dist = linspace(mindist,maxdist,nbins_dist+1);
centers_dist = (edges_dist(1:end-1)+edges_dist(2:end))/2;

fracperfly = nan([nbins_dist,nflies]);
for i = 1:nflies,
  % rows are distance, columns are angle
  counts = hist(dcenter{i},centers_dist);
  fracperfly(:,i) = counts / sum(counts(:));
end

if ishandle(1),
  set(0,'CurrentFigure',1);
else
  figure(1);
end
clf;
idx = wingtypeidx == 1;
hfly = plot(centers_dist,fracperfly(:,idx),'.-','Color',[.7,.7,.7]); %#ok<NASGU>
hold on;
meanfrac = mean(fracperfly(:,idx),2);
hmean = plot(centers_dist,meanfrac,'k.-','LineWidth',2); %#ok<NASGU>

axisalmosttight;
drawnow;
