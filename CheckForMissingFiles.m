function [ismissingrequiredfile,missingrequiredfiles,...
  ismissingdesiredfile,missingdesiredfiles] = ...
  CheckForMissingFiles(expdir,cbparams,stage)

ismissingrequiredfile = false;
ismissingdesiredfile = false;
missingrequiredfiles = {};
missingdesiredfiles = {};

files = fieldnames(cbparams.dataloc);
isrequiredstage = structfun(@(x) isstruct(x) && isfield(x,'essential') && ...
  (x.essential>=1) && isfield(x,'type') && strcmp(x.type,stage),cbparams.dataloc);
required_files = files(isrequiredstage);
isdesiredstage = structfun(@(x) isstruct(x) && isfield(x,'essential') && ...
  (x.essential==0) && isfield(x,'type') && strcmp(x.type,stage),cbparams.dataloc);
desired_files = files(isdesiredstage);

for i = 1:numel(required_files),
  file = cbparams.dataloc.(required_files{i});
  if isfield(file,'searchstr'),
    filestr = file.searchstr;
  else
    filestr = file.filestr;
  end 
  if any(filestr == '*'),
    isfile = ~isempty(dir(fullfile(expdir,filestr)));
  else
    isfile = exist(fullfile(expdir,filestr),'file');
  end
  
  if ~isfile,
    missingrequiredfiles{end+1} = required_files{i}; %#ok<AGROW>
    ismissingrequiredfile = true;
  end
end

if strcmpi(stage,'compute_perframe_features'),
  perframefns = cbparams.compute_perframe_features.perframefns;
    
  perframedir = fullfile(expdir,cbparams.dataloc.perframedir.filestr);
  for i = 1:numel(perframefns),
    if ~exist(fullfile(perframedir,[perframefns{i},'.mat']),'file'),
      missingrequiredfiles{end+1} = ['perframe_',perframefns{i}]; %#ok<AGROW>
    end
  end  
  
end


for i = 1:numel(desired_files),
  file = cbparams.dataloc.(desired_files{i});
  if isfield(file,'searchstr'),
    filestr = file.searchstr;
  else
    filestr = file.filestr;
  end 
  if any(filestr == '*'),
    isfile = ~isempty(dir(fullfile(expdir,filestr)));
  else
    isfile = exist(fullfile(expdir,filestr),'file');
  end
  
  if ~isfile,
    missingdesiredfiles{end+1} = desired_files{i}; %#ok<AGROW>
    ismissingdesiredfile = true;
  end
end