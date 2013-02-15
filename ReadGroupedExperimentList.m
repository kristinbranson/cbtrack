function expdirs = ReadGroupedExperimentList(expfile)

fid = fopen(expfile,'r');
endread = false;
while true,
  
  isfirst = true;
  while true,
    l = fgetl(fid);
    if ~ischar(l),
      endread = true;
      break;
    end
    l = strtrim(l);
    if isempty(l),
      break;
    end
    if isfirst,
      conditionnamecurr = l;
      expdirscurr = {};
      isfirst = false;      
    else
      expdirscurr{end+1} = l;
    end    
  end
  
  if ~isfirst,
    expdirs.(conditionnamecurr) = expdirscurr;
  end
  
  if endread,
    break;
  end
end

fclose(fid);