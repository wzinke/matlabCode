function hhh=hline(y,varargin)
% function h=hline(y, linetype, label)
% 
% Modified by Kaleb Lowe 6/25/15
%    -Changed default linetype to '-k'. 
%    -Cut out multi-line support
%    -Changed inputs to varargin, to set output parameters according to
%     normal line specifications
%
% Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
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
% h = hline(42,'g','The Answer')
%
% returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
% the current axes, close to the line, which reads "The Answer".
%
% hline also supports vector inputs to draw multiple lines at once.  For example,
%
% hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
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

if size(varargin,1) == length(y)
    uniform = 0;
else
    uniform = 1;
end

x=get(gca,'xlim');
for iy = 1:length(y)
    h=plot(x,[y(iy) y(iy)],linetype);
    for iv = 1:2:size(varargin,2)
        if uniform
            if ~isempty(varargin{1,iv}),
                set(h,varargin{1,iv},varargin{1,iv+1});
            end
        else
            if ~isempty(varargin{iy,iv}),
                set(h,varargin{iy,iv},varargin{iy,iv+1});
            end
        end
    end
end

if g==0
hold off
end
set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends

if nargout
    hhh=h;
end
