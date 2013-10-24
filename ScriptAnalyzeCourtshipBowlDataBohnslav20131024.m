%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/simplewing;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/compute_perframe_features;

%% parameters

expdir = '/groups/branson/bransonlab/projects/CourtshipBowls/data/Bohnslav/WT6dpf30deg_2013-10-16-152942-0000_arena_3';
moviefilestr = 'movie.avi';
dofixbg = false;
analysis_protocol = '20131024_bohnslav';
settingsdir = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/settings';
paramsfilestr = 'params.xml';
cbparams = ReadXMLParams(fullfile(fullfile(settingsdir,analysis_protocol),paramsfilestr));

params = cbparams.track;

SetBackgroundTypes;
if ischar(params.bgmode) && isfield(bgtypes,params.bgmode),
  params.bgmode = bgtypes.(params.bgmode);
end
params.DEBUG = 2;
if params.dotrackwings || strcmp(params.assignidsby,'wingsize'),
  params.wingtracking_params = cbparams.wingtrack;
end

%% open video

moviefile = fullfile(expdir,moviefilestr);
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
im = readframe(1);
[nr,nc,ncolors] = size(im);

%% create dummy roi data

roidata = struct;
roidata.nrois = 1;
roidata.roibbs = [1,nc,1,nr];
roidata.inrois = {true(nr,nc)};
roidata.nflies_per_roi = 2;


%% estimate bg model

[bgmed,bgfixdata] = EstimateBGModel(readframe,nframes,cbparams.track,dofixbg);


%% set bgsub thresholds
% 
% done = false;
% f = 1;
% fold = nan;
% 
% 
% hfig = 1;
% figure(hfig);
% clf;
% hax = gca;
% 
% expi = 0;
% fid = 0;
% 
% him = nan;
% h1 = nan;
% h2 = nan;
% hti = nan;
% 
% bgthresh = cbparams.track.bgthresh;
% bgthresh_low = cbparams.track.bgthresh_low;
% 
% while true,
% 
%   [~,experiment_name] = fileparts(expdir);
%   imheight = headerinfo.nr;
%   imwidth = headerinfo.nc;
%     
%   if f ~= fold,
%     
%     im = readframe(f);
%     if size(im,3) > 1,
%       im = rgb2gray(im);
%     end
%     
%     % subtract off background
%     switch cbparams.track.bgmode,
%       case 'DARKBKGD',
%         dbkgd = imsubtract(im,bgmed);
%       case 'LIGHTBKGD',
%         dbkgd = imsubtract(bgmed,im);
%       case 'OTHERBKGD',
%         dbkgd = imabsdiff(im,bgmed);
%     end
%     
%     fold = f;
% 
%   end
%   
%   % threshold
%   isfore = dbkgd >= bgthresh;
%   b = bwboundaries(isfore);
%   rc = zeros(0,2);
%   for i = 1:numel(b),
%     rc = [rc;nan(1,2);b{i}];
%   end
% 
%   isforelow = dbkgd >= bgthresh_low;
%   b = bwboundaries(isforelow);
%   rclow = zeros(0,2);
%   for i = 1:numel(b),
%     rclow = [rclow;nan(1,2);b{i}];
%   end
% 
%   
%   if ishandle(him),
%     set(him,'CData',im);
%   else
%     hold(hax,'off');
%     him = imagesc(im,'Parent',hax,[0,255]);
%     axis(hax,'image');
%     colorbar('Peer',hax);
%   end
%   if ishandle(h1),
%     set(h1,'XData',rclow(:,2),'YData',rclow(:,1));
%   else
%     hold(hax,'on');
%     h1 = plot(hax,rclow(:,2),rclow(:,1),'b.');
%   end
%   if ishandle(h2),
%     set(h2,'XData',rc(:,2),'YData',rc(:,1));
%   else
%     hold(hax,'on');
%     h2 = plot(hax,rc(:,2),rc(:,1),'r.');
%   end
%   if ishandle(hti),
%     set(hti,'String',num2str(f));
%   else
%     hti = title(hax,num2str(f));
%   end
%   
%   
%   while true,
%     
%     fprintf('1: decrease low threshold\n');
%     fprintf('2: increase low threshold\n');
%     fprintf('3: decrease high threshold\n');
%     fprintf('4: increase high threshold\n');
%     fprintf('5: previous frame\n');
%     fprintf('6: next frame\n');
%     fprintf('7: select frame\n');
%     fprintf('8: select random frame\n');
%     fprintf('0: exit\n');
%     
%     v = input('');
%     
%     switch v,
%       case 1,
%         bgthresh_low = bgthresh_low*.99;
%       case 2,
%         bgthresh_low = bgthresh_low/.99;
%       case 3,
%         bgthresh = bgthresh*.99;
%       case 4,
%         bgthresh = bgthresh/.99;
%       case 5,
%         f = max(1,f-1);
%       case 6,
%         f = min(nframes,f+1);
%       case 7,
%         while true,
%           f1 = input(sprintf('Enter frame (1 to %d)',nframes));
%           if f1 < 1 || f1 > nframes || round(f1) ~= f1,
%             continue;
%           end
%           break;
%         end
%         f = f1;
%       case 8,
%         f = randsample(nframes,1);        
%       case 0,
%         done = true;
%         break;
%       otherwise,
%         continue;
%     end
%     
%     break;
%     
%   end
% 
%   if done,
%     break;
%   end
%   
% end



%% track

tmpfilename = fullfile(expdir,'tmptracking.mat');
trackdata = TrackTwoFlies(moviefile,bgmed,roidata,params,'tmpfilename',tmpfilename,'logfid',1);

save(fullfile(expdir,cbparams.dataloc.trx.filestr),'-struct','trackdata');

%% results movie

CourtshipBowlMakeResultsMovie(expdir,'analysis_protocol',analysis_protocol);
