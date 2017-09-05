function hhh=vline(x,varargin)
% function h=vline(x, linetype, label)
% 
% Modified by Kaleb Lowe 06/25/15
%    - Changed default linetype to '-k'
%    - Cut out multi-line support
%    - 'x' is the only required input now, varargin can specify line
%      properties similar to usual plotting schemes
%    - Re-implemented multi-line support in a different way

% Draws a vertical line on the current axes at the location specified by 'x'.  Optional arguments are
% 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
% label appears in the same color as the line.
%
% The line is held on the current axes, and after plotting the line, the function returns the axes to
% its prior hold state.
%
% The HandleVisibility property of the line object is set to "off", so not only does it not appear on
% legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
% return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
% overridden by setting the root's ShowHiddenHandles property to on.
%
% h = vline(42,'g','The Answer')
%
% returns a handle to a green vertical line on the current axes at x=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% vline also supports vector inputs to draw multiple lines at once.  For example,
%
% vline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
%
% draws three lines with the appropriate labels and colors.
% 
% By Brandon Kuczenski for Kensington Labs.
% brandon_kuczenski@kensingtonlabs.com
% 8 November 2001

% Set defaults
linetype='-k';
label='';
   
g=ishold(gca);
hold on

if unique(size(varargin)) == 1,
    varargin = varargin{1};
else
    uniform = 1;
end

if size(varargin,1) == length(x)
    uniform = 0;
else
    uniform = 1;
end

y=get(gca,'ylim');
%y = [-inf inf];
for ix = 1:length(x)
    h=plot([x(ix) x(ix)],y,linetype);
    for iv = 1:2:length(varargin)
        if uniform
            if ~isempty(varargin{1,iv})
                set(h,varargin{1,iv},varargin{1,iv+1});
            end
        else
            if ~isempty(varargin{ix,iv}),
                set(h,varargin{ix,iv},varargin{ix,iv+1});
            end
        end
    end
end

if g==0
hold off
end
set(h,'tag','vline','handlevisibility','off')

if nargout
    hhh=h;
end
