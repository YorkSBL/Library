function xy = extractall(ah,properties)
% function xy = extractall(?axis_or_figure_handle=gca?,?properties=false?)
%  
% Extract all xy data (lines) from a figure. Returns a structure
% array containing data values (.x,.y) and optional property values 
% (.info) for the data sets displayed in the figure with the given axis 
% (or figure) handle (gca by default). The nth data set is in xy(n).  
% extractall() is useful for extracting data from figures.
% CASHera  

  if (nargin < 1 | ~ishandle(ah) | ~strcmp (get(ah,'Type'),'axes'))
    ah = gca;
  end
  if (nargin < 2)
    properties = false;
  end
  
  xy = [];

  % sort by handle to put them in the order in which they were added...
  lhs = sort(findall(ah,'Type','line'));
  
  N = numel(lhs);
  for n=1:N
    xy(n).x = get(lhs(n),'XData');
    xy(n).y = get(lhs(n),'YData');
    if (properties)
      xy(n).info = get(lhs(n));
    end
  end
  
  return