function filename = readunixlinks(filename)

if ~isunix,
  return;
end
[~,link] = unix(sprintf('readlink %s',filename));
link = strtrim(link);
if isempty(link),
  return;
end

% is it relative?
if link(1) ~= filesep,
  [pathcurr] = fileparts(filename);
  link = fullfile(pathcurr,link);
end

filename = readunixlinks(link);
