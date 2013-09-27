function [bgmed,bgfixdata] = FixBgModel(bgmed,moviefile,tracking_params)

bgmed0 = bgmed;
hfig = figure;
hax = gca;
him = imagesc(bgmed,'Parent',hax,[0,255]);
colormap(hax,'gray');
axis(hax,'image');
hold(hax,'on');
hroi = [];

hplayfmf = playfmf('moviefile',moviefile);
options = struct;
options.Resize='on';

startframe = 1;
endframe = 1;

[readframe,nframes,fid] = get_readframe_fcn(moviefile);

[nr,nc,ncolors] = size(bgmed);
if ncolors > 1,
  error('Not implemented for color images yet');
end
[XGRID,YGRID] = meshgrid(1:nc,1:nr);

bgfixdata = struct;
bgfixdata.xs = {};
bgfixdata.ys = {};

res = questdlg('Does the background model need fixing?');
if ~strcmp(res,'Yes'),
  return;
end

while true,
  
  bgmedprev = bgmed;
  title(hax,'Click to select a region of interest. Press ENTER when done.');
  
  [x,y] = getline(hfig,'closed');
  hroi(end+1) = plot(hax,x,y,'.-','Color',[.8,0,0],'LineWidth',2);
  
  if ~ishandle(hplayfmf),
    hplayfmf = playfmf('moviefile',moviefile);
  end
  handles = guidata(hplayfmf);

  holdstate = ishold(handles.axes_Video);

  hold(handles.axes_Video,'on');
  hextra = plot(handles.axes_Video,x,y,'.-','Color',[.8,0,0],'LineWidth',2);
  if ~holdstate,
    hold(handles.axes_Video,'off');
  end

  while true,
    startframe = input(sprintf('Start of interval for which there is no fly in this ROI (number between 1 and %d): ',nframes));
    if isempty(startframe) || round(startframe) ~= startframe || startframe < 1 || startframe > nframes,
      fprintf('Start frame <= end frame must be integers between 1 and %d',nframes);
      continue;
    end
    break;
  end
  while true,
    endframe = input(sprintf('End of interval for which there is no fly in this ROI (number between %d and %d): ',startframe,nframes));
    if isempty(endframe) || round(endframe) ~= endframe || endframe < startframe || endframe > nframes,
      fprintf('Start frame <= end frame must be integers between 1 and %d',nframes);
      continue;
    end
    break;
  end
  
  inroi = inpolygon(XGRID,YGRID,x,y);
  if endframe - startframe + 1 <= tracking_params.bg_nframes,
    fs = startframe:endframe;
  else
    fs = unique(round(linspace(startframe,endframe,tracking_params.bg_nframes)));
  end
  
  buffer = readframe(1);
  buffer = buffer(inroi);
  buffer = repmat(buffer,[1,numel(fs)]);
  for i = 1:numel(fs),
    f = fs(i);
    tmp = readframe(f);
    buffer(:,i) = tmp(inroi);
  end
  bgmedcurr = uint8(median(single(buffer),2));

  bgmed(inroi) = bgmedcurr;
  
  if ishandle(hextra),
    delete(hextra);
  end
  
  bgfixdata.xs{end+1} = x;
  bgfixdata.ys{end+1} = y;
  
  set(him,'CData',bgmed);
  
  res = questdlg('What next?','What next?','Add a new ROI','Undo','Finished','Add a new ROI');
  if strcmp(res,'Finished'),
    break;
  elseif strcmp(res,'Undo'),
    bgmed = bgmedprev;
    if ishandle(hroi(end)),
      delete(hroi(end));
      hroi(end) = [];
    end
    bgfixdata.xs(end) = [];
    bgfixdata.ys(end) = [];
    set(him,'CData',bgmed);
  end    
  
end

if fid > 0,
  fclose(fid);
end
