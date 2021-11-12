function [y,x,yl,yu] = extract(ah,viewflag)
% function [y,x,yl,yu] = extract(?axis_or_figure_handle=gca?,?viewflag='x'?)
%  
% Extract a data selection from a figure. Returns coordinate 
% values (y,x) for the data set currently displayed in the 
% the figure having the given axis (or figure) handle (gca by default).
% Extracts only those data points shown in the current 
% (possibly zoomed) view according to the optional viewflag, 
% which  indicates whether to define the data selection by 
% the current x-axis limis ('x'), the y-axis limits ('y'), 
% or by both ('xy'). The default is 'x'.  
%  
% If more than one data set appears in the view, extract() returns 
% data for the line currently selected (e.g., using the mouse).  
% [If more than one line is selected---or if no lines
% are selected at all---extract() uses the first 
% such line encountered in the list of axis children.]
%
% Extract() is useful for extracting data from figures and for
% making data selections after zooming.  
% CASHera  

  if (nargin<1 | ~ishandle(ah) | ~strcmp (get(ah,'Type'),'axes'))
    ah = gca;
  end
  if (nargin<2)
    viewflag = 'x';
  end
  
  % find all lines and note whether they're selected...
  nl = 0;				% number of lines found
  selected = [];
  Children = get(ah,'Children'); 
  for n=1:length(Children)
    child = Children(n);
    if (strcmp (get(child,'Type'),'line') | strcmp (get(child,'Type'),'hggroup'))
      nl = nl + 1;
      Lines(nl).x = get(child,'XData');
      Lines(nl).y = get(child,'YData');
      try
        Lines(nl).yl = get(child,'LData');
        Lines(nl).yu = get(child,'UData');
      catch
        Lines(nl).yl = [];
        Lines(nl).yu = [];
      end
      selected(nl) = strcmp(get(child,'Selected'),'on');
    end
  end

  % make sure we've got at least one line...
  if (nl == 0)
    x = []; y = []; yl = []; yu = [];
    warning('No lines in figure.');
    return
  end

  % use only selected lines (or all if none selected)...
  idxs = eitheror(any(selected),find(selected),1:nl);

  % find first line with data in view...
  xlim = get(ah,'XLim');
  ylim = get(ah,'YLim');
  for idx=idxs
    x = Lines(idx).x;
    y = Lines(idx).y;

    okx = ininterval(x,xlim,2);
    oky = ininterval(y,ylim,2);
  
    switch (viewflag)
     case {'x'}
      ok = okx;
     case {'y'}
      ok = oky;
     case {'xy','yx'}
      ok = oky & oky;
     otherwise
      ok = okx;
    end
    
    if (any(ok))
      break
    end
  end

  if (~any(ok))
    x = []; y = []; yl = []; yu = [];
    warning ('No data in view.');
    return
  end

  x = x(ok);
  y = y(ok);
  
  if (~isempty(Lines(idx).yl))
    yl = Lines(idx).yl(ok);
    yu = Lines(idx).yu(ok);
  end
  
  return