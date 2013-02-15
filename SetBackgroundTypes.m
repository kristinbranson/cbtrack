% SetBackSubTypes

DARKBKGD = 1;
LIGHTBKGD = -1;
OTHERBKGD = 0;

bgtypes.DARKBKGD = DARKBKGD;
bgtypes.LIGHTBKGD = LIGHTBKGD;
bgtypes.OTHERBKGD = OTHERBKGD;

% if isfield(tracking_params,'bgmode'),
%   SetBackgroundTypes;
%   if ischar(tracking_params.bgmode) && isfield(bgtypes,tracking_params.bgmode),
%     tracking_params.bgmode = bgtypes.(tracking_params.bgmode);
%   end
% end