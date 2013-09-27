function [meannormfrac,meanfrac,normfracperfly,fracperfly,edges_angle,edges_dist] = ...
  HistogramRelativePositions(expdir,varargin)

nbins_angle = 12;
nbins_dist = 10;
mindist = 1;
maxdist = 27;

[nbins_angle,nbins_dist,mindist,maxdist] = myparse(varargin,...
  'nbins_angle',nbins_angle,...
  'nbins_dist',nbins_dist,...
  'mindist',mindist,...
  'maxdist',maxdist);

load(fullfile(expdir,'registered_trx.mat'),'trx');
nflies = numel(trx);
wingtype = cell(1,nflies);
for i = 1:nflies,
  wingtype{i} = trx(i).wingtype{1};
end

[unique_wingtypes,~,wingtypeidx] = unique(wingtype);

load(fullfile(expdir,'perframe','dcenter.mat'),'data','units');
dcenter = data;
load(fullfile(expdir,'perframe','anglefrom1to2_nose2ell.mat'),'data','units');
anglefrom1to2 = data;
for i = 1:nflies,
  anglefrom1to2{i} = anglefrom1to2{i};
end

edges_angle = linspace(-pi,pi,nbins_angle+1);
centers_angle = (edges_angle(1:end-1)+edges_angle(2:end))/2;

%edges_dist = logspace(log10(mindist),log10(maxdist),nbins_dist+1);
edges_dist = linspace(mindist,maxdist,nbins_dist+1);
centers_dist = (edges_dist(1:end-1)+edges_dist(2:end))/2;

fracperfly = nan([nbins_dist,nbins_angle,nflies]);
for i = 1:nflies,
  % rows are distance, columns are angle
  counts = hist3([dcenter{i}(:),anglefrom1to2{i}(:)],{centers_dist,centers_angle});
  %counts(end,1) = max(counts(:)*2);
  fracperfly(:,:,i) = counts / sum(counts(:));
end

binarea = pi*diff(edges_dist.^2)/2/nbins_angle;
normfracperfly = bsxfun(@rdivide,fracperfly,binarea(:));

figure(1);
nr = 2;
nc = nflies/2;
clf;
hax = createsubplots(nr,nc,[0,.01]);
h = nan(1,nflies);
for i = 1:nflies,
  axes(hax(i));
  h(i) = polarimagesc(edges_dist,edges_angle+pi/2,[0,0],normfracperfly(:,:,i));
  title(wingtype{i});
  axis equal off;
  axisalmosttight;
  %set(hax(i),'CLim',[min(normfracperfly(:)),max(normfracperfly(:))]);
end

figure(2);
clf;
hax = createsubplots(1,numel(unique_wingtypes),[0,.01]);
meanfrac = struct;
meannormfrac = struct;
for i = 1:numel(unique_wingtypes),
  axes(hax(i));
  wingtypecurr = unique_wingtypes{i};
  meanfrac.(wingtypecurr) = mean(fracperfly(:,:,wingtypeidx==i),3);
  meannormfrac.(wingtypecurr) = mean(normfracperfly(:,:,wingtypeidx==i),3);
  polarimagesc(edges_dist,edges_angle+pi/2,[0,0],meannormfrac.(wingtypecurr));
  title(wingtypecurr);
  axis equal off;
  axisalmosttight;
end
drawnow;
