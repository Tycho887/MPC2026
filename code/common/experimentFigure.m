function fig = experimentFigure(varargin)
%EXPERIMENTFIGURE Create a light figure for legible report graphics.
% Set the theme before creating axes, without changing the user's defaults.
fig = figure(varargin{:});
if isprop(fig,'Theme'), set(fig,'Theme','light'); end
set(fig,'Color','w');
end
