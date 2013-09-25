rootdatadir = '/groups/branson/bransonlab/projects/FlyBowl/Courtship_plates_movies';
experiment_names = {
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlA_20130207T135248'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlA_20130207T145947'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlB_20130207T135253'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlB_20130207T145954'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlC_20130207T135756'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlC_20130207T150512'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlD_20130207T135802'
  'EXT_CantonS_1220002_TrpA_Rig2Plate17BowlD_20130207T150518'
  };
handling_protocol_dict = {
  'HP_Olympiad_v007p2.xls','rejected'
  'HP_Olympiad_v007p1.xls','mated'
  'HP_Olympiad_v007p3.xls','naive'
  'HP_Olympiad_v007p7.xls',{'mated_normal_rejected_wingclipped','mated_wingclipped_rejected_normal'}};

conditions = cell(1,numel(experiment_names));
for i = 1:numel(experiment_names),
  expdir = fullfile(rootdatadir,experiment_names{i});
  metadata = ReadMetadataFile(fullfile(expdir,'Metadata.xml'));
  j = find(strcmp(handling_protocol_dict(:,1),metadata.handling_protocol));
  if isempty(j),
    error('Could not find protocol %s',metadata.handling_protocol);
  end
  if ischar(handling_protocol_dict{j,2}),
    conditions{i} = handling_protocol_dict{j,2};
  else
    while true,
      fprintf('%s\n',experiment_names{i});
      fprintf(['Notes technical: ',metadata.notes_technical,'\n']);
      for k = 1:numel(handling_protocol_dict{j,2}),
        fprintf('Enter %d for %s\n',k,handling_protocol_dict{j,2}{k});
      end
      res = input('');
      if res > 0 && res <= numel(handling_protocol_dict{j,2}) && round(res)==res,
        k = res;
        break;
      end
    end
    conditions{i} = handling_protocol_dict{j,2}{k};
    fprintf('Chose %s\n\n',handling_protocol_dict{j,2}{k});
  end
end

fid = fopen('expdirlists/grouped_expdirs_galit_rejection_wingclipped_20130207.txt','w');
[unique_conditions,~,idx] = unique(conditions);
for i = 1:numel(unique_conditions),
  fprintf(fid,'%s\n',unique_conditions{i});
  for j = find(idx==i),
    fprintf('%s\n',fullfile(rootdatadir,experiment_names{j}));
  end
  fprintf('\n');
end
fclose(fid);