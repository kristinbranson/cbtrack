% make results movies
function resultsmoviedata = CourtshipBowlMakeResultsMovie(expdir,varargin)

version = '0.1.1';
timestamp = datestr(now,TimestampFormat);

resultsmoviedata = struct;
resultsmoviedata.cbresultsmovie_version = version;
resultsmoviedata.cbresultsmovie_timestamp = timestamp;

ParseCourtshipBowlParams;
[DEBUG] = myparse(leftovers,'debug',false);



%% movie parameters

resultsmovie_params = cbparams.results_movie;
defaulttempdatadir = '/groups/branson/bransonlab/projects/CourtshipBowls/TempData_CourtshipBowlMakeResultsMovie';
if ~isfield(resultsmovie_params,'tempdatadir'),
  resultsmovie_params.tempdatadir = defaulttempdatadir;
elseif isunix
  [status1,res] = unix(sprintf('echo %s',resultsmovie_params.tempdatadir));
  if status1 == 0,
    resultsmovie_params.tempdatadir = strtrim(res);
  end
end

if ~exist(resultsmovie_params.tempdatadir,'dir'),
  [success1,msg1] = mkdir(resultsmovie_params.tempdatadir);
  if ~success1,
    error('Error making directory %s: %s',resultsmovie_params.tempdatadir,msg1);
  end
end

%% log file

if isfield(cbparams.dataloc,'resultsmovie_log') && ~DEBUG,
  logfile = fullfile(expdir,cbparams.dataloc.resultsmovie_log.filestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

fprintf(logfid,'\n\n***\nRunning CourtshipBowlMakeResultsMovie version %s analysis_protocol %s at %s\n',version,real_analysis_protocol,timestamp);
resultsmoviedata.analysis_protocol = real_analysis_protocol;

%% location of data

[~,basename] = fileparts(expdir);
moviefile = fullfile(expdir,cbparams.dataloc.movie.filestr);
trxfile = fullfile(expdir,cbparams.dataloc.trx.filestr);
avifilestr = sprintf('%s_%s',cbparams.dataloc.resultsavi.filestr,basename);
avifile = fullfile(resultsmovie_params.tempdatadir,[avifilestr,'_temp.avi']);
xvidfile = fullfile(expdir,[avifilestr,'.avi']);

%% read start and end of cropped trajectories

fprintf(logfid,'Reading in trajectories...\n');

load(trxfile,'trx');
end_frame = max([trx.endframe]);
start_frame = min([trx.firstframe]);
nframes = end_frame-start_frame + 1;
nflies = numel(trx);
firstframes_off = min(max(0,round(resultsmovie_params.firstframes*nframes)),nframes-1);
firstframes_off(resultsmovie_params.firstframes < 0) = nan;
middleframes_off = round(resultsmovie_params.middleframes*nframes);
middleframes_off(resultsmovie_params.middleframes < 0) = nan;
endframes_off = round(resultsmovie_params.endframes*nframes);
endframes_off(resultsmovie_params.endframes < 0) = nan;
idx = ~isnan(middleframes_off);
firstframes_off(idx) = ...
  min(nframes-1,max(0,middleframes_off(idx) - ceil(resultsmovie_params.nframes(idx)/2)));
idx = ~isnan(endframes_off);
firstframes_off(idx) = ...
  min(nframes-1,max(0,endframes_off(idx) - resultsmovie_params.nframes(idx)));
endframes_off = firstframes_off + resultsmovie_params.nframes - 1;


firstframes = start_frame + firstframes_off;

%% option to not specify nzoomr, nzoomc

fprintf(logfid,'Setting parameters...\n');

if ischar(resultsmovie_params.nzoomr) || ischar(resultsmovie_params.nzoomc),
    
  if isnumeric(resultsmovie_params.nzoomr),
    nzoomr = resultsmovie_params.nzoomr;
    nzoomc = round(nflies/nzoomr);
  elseif isnumeric(resultsmovie_params.nzoomc),
    nzoomc = resultsmovie_params.nzoomc;
    nzoomr = round(nflies/nzoomc);
  else
    nzoomr = ceil(sqrt(nflies));
    nzoomc = round(nflies/nzoomr);
  end
  resultsmovie_params.nzoomr = nzoomr;
  resultsmovie_params.nzoomc = nzoomc;
  
  if iscell(resultsmovie_params.figpos),  
    [readframe,~,fid] = get_readframe_fcn(moviefile);
    im = readframe(1);
    [nr,nc,~] = size(im);
    
    rowszoom = floor(nr/nzoomr);
    imsize = [nr,nc+rowszoom*nzoomc];
    figpos = str2double(resultsmovie_params.figpos);
    if isnan(figpos(3)),
      figpos(3) = figpos(4)*imsize(2)/imsize(1);
    elseif isnan(figpos(4)),
      figpos(4) = figpos(3)*imsize(1)/imsize(2);
    end
    resultsmovie_params.figpos = figpos;
    
    if fid > 1,
      fclose(fid);
    end
  end
  
  
  
end

%% choose colors for each fly

fprintf(logfid,'Choosing colors for each fly...\n');

colors = jet(nflies)*.8;

% alone flies go in the middle
isalone = false(1,nflies);
maxroi = max([trx.roi]);
for roii = 1:maxroi,
  idx = [trx.roi] == roii;
  if nnz(idx) == 1,
    isalone(idx) = true;
  end
end

nalone = nnz(isalone);
npairs = (nflies - nalone)/2;

colorsp1 = colors(1:npairs,:);
colorsalone = colors(npairs+1:npairs+nalone,:);
colorsp2 = colors(npairs+nalone+1:end,:);

ip = 1;
ialone = 1;
for roii = 1:maxroi,
  idx = find([trx.roi] == roii);
  if nnz(idx) == 0,
    continue;
  elseif numel(idx) == 1,
    colors(idx,:) = colorsalone(ialone,:);
    ialone = ialone+1;
  else
    colors(idx(1),:) = colorsp1(ip,:);
    colors(idx(2),:) = colorsp2(ip,:);
    ip = ip+1;
  end
end

%% create movie

fprintf(logfid,'Calling make_ctrax_results_movie...\n');

if ~DEBUG && exist(avifile,'file'),
  try
    delete(avifile);
  catch ME,
    fprintf(logfid,'Could not remove avi file %s:\n%s\n',avifile,getReport(ME));
  end
end

[succeeded,~,~,height,width]= ...
  make_ctrax_result_movie('moviename',moviefile,'trxname',trxfile,'aviname',avifile,...
  'nzoomr',resultsmovie_params.nzoomr,'nzoomc',resultsmovie_params.nzoomc,...
  'boxradius',resultsmovie_params.boxradius,'taillength',resultsmovie_params.taillength,...
  'fps',resultsmovie_params.fps,...
  'maxnframes',resultsmovie_params.nframes,...
  'firstframes',firstframes,...
  'figpos',resultsmovie_params.figpos,...
  'movietitle',basename,...
  'compression','none',...
  'useVideoWriter',false,...
  'titletext',false,...
  'avifileTempDataFile',[avifile,'-temp'],...
  'dynamicflyselection',true,...
  'doshowsex',true,...
  'colors',colors);

if ishandle(1),
  close(1);
end

if ~succeeded,
  error('Failed to create raw avi %s',avifile);
end

%% create subtitle file

fprintf(logfid,'Creating subtitle file...\n');

subtitlefile = fullfile(expdir,'subtitles.srt');
fid = fopen(subtitlefile,'w');
dt = [0,resultsmovie_params.nframes];
ts = cumsum(dt);
for i = 1:numel(dt)-1,
  fprintf(fid,'%d\n',i);
  fprintf(fid,'%s --> %s\n',...
    datestr(ts(i)/resultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'),...
    datestr((ts(i+1)-1)/resultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'));
  fprintf(fid,'%s, fr %d-%d\n\n',basename,...
    firstframes_off(i)+1,...
    endframes_off(i)+1);
end
fclose(fid);

%% compress

fprintf(logfid,'Compressing to xvid avi file...\n');

tmpfile = [xvidfile,'.tmp'];
newheight = 4*ceil(height/4);
newwidth = 4*ceil(width/4);
% subtitles are upside down, so encode with subtitles and flip, then flip
% again
cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf scale=%d:%d,flip -sub %s -subfont-text-scale 2 -msglevel all=2',...
  avifile,tmpfile,newwidth,newheight,subtitlefile);
status = system(cmd);
if status ~= 0,
  fprintf('*****\n');
  warning('Failed to compress avi to %s',xvidfile);
  fprintf('Need to run:\n');
  fprintf('%s\n',cmd);
  cmd2 = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
    tmpfile,xvidfile);
  fprintf('then\n');
  fprintf('%s\n',cmd2);
  fprintf('then delete %s %s %s\n',tmpfile,avifile,subtitlefile);
  fprintf('*****\n');
else
  cmd = sprintf('mencoder %s -o %s -ovc xvid -xvidencopts fixed_quant=4 -vf flip -msglevel all=2',...
    tmpfile,xvidfile);
  status = system(cmd);
  if status ~= 0,
    fprintf('*****\n');
    warning('Failed to add subtitles to %s',xvidfile);
    fprintf('Need to run:\n');
    fprintf('%s\n',cmd);
    fprintf('then delete %s %s %s\n',tmpfile,avifile,subtitlefile);
    fprintf('*****\n');    
  else
    delete(tmpfile);
    delete(avifile);
    delete(subtitlefile);
  end
end

%% save info to mat file

resultsmoviedata.resultsmovie_params = resultsmovie_params;
resultsmoviedata.firstframes = firstframes;

filename = fullfile(expdir,cbparams.dataloc.resultsmoviedatamat.filestr);
fprintf(logfid,'Saving info to mat file %s...\n',filename);

if exist(filename,'file'),
  delete(filename);
end
save(filename,'-struct','resultsmoviedata');

%% close log

fprintf(logfid,'Finished running CourtshipBowlMakeResultsMovie at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end
