function deplot_callbacks (action)
% GUI callbacks for deplot.m
% C.A. Shera
  global DEPLOT;

  switch (action)
   case 'quit'				% quit deplot...
    close (gcbf);
   case 'help'				% popup some help...
    helpdlg ({'Keyboard commands when digitizing:', ...
	      '    z = zoom in', ...
	      '    u = unzoom (zoom out)', ...
	      '    U = unzoom all the way', ...
	      '    n = begin new curve', ...
	      '    s = show current coordinates', ...
	      '    q = quit digitizing'},'dePlot Help');
   case 'load'				% load a scanned image...
    if (load_image)
      set (findobj(gcbf,'Tag','DigitizeButton'),'Enable','on');
    end    
   case 'save'				% save digitized data...
    save_data;
   case 'digitize'			% digitize...
    enable_buttons ('off');
    if (~read_values)
      enable_buttons ('on');
      return
    end
    if (~digitize_axes)
      % erase the title and re-enable the buttons...
      title ('','Color',[1,0,0]);
      enable_buttons ('on');
      return
    end
    if (~find_matrix)
      enable_buttons ('on');
      return
    end
    DEPLOT.data = digitize_points;
    % erase the title and re-enable the buttons...
    title ('','Color',[1,0,0]);
    enable_buttons ('on');
end

function status = read_values
  global DEPLOT;
  status = 1;
  
  h = findobj(gcbf,'Tag','XLogCheckbox');
  DEPLOT.x_log = get (h,'Value');

  h = findobj(gcbf,'Tag','YLogCheckbox');
  DEPLOT.y_log = get (h,'Value');

  h = findobj(gcbf,'Tag','OriginCheckbox');
  DEPLOT.explicit_origin = get(h,'Value');

  h = findobj(gcbf,'Tag','PerpAxesCheckbox');
  DEPLOT.force_perp_axes = get (h,'Value');

  h = findobj(gcbf,'Tag','DrawAxesCheckbox');
  DEPLOT.draw_axes = get (h,'Value');

  eh = findobj(gcbf,'Tag','XminEditText');
  x_min = str2num(get(eh,'String'));

  eh = findobj(gcbf,'Tag','XmaxEditText');
  x_max = str2num(get(eh,'String'));

  eh = findobj(gcbf,'Tag','YminEditText');
  y_min = str2num(get(eh,'String'));

  eh = findobj(gcbf,'Tag','YmaxEditText');
  y_max = str2num(get(eh,'String'));

  % perform a few basic sanity checks...
  if (DEPLOT.x_log)
    if (x_min <= 0 | x_max <= 0)
      uiwait (errordlg ('Illegal limits for logarithmic x axis.'));
      status = 0;
      return
    else
      x_min = log (x_min);
      x_max = log (x_max);
    end
  end
  if (DEPLOT.y_log)
    if (y_min <= 0 | y_max <= 0)
      uiwait (errordlg ('Illegal limits for logarithmic y axis.'));
      status = 0;
      return
    else
      y_min = log (y_min);
      y_max = log (y_max);
    end
  end
  
  if (abs(x_max-x_min)<1e-10 | abs(y_max-y_min)<1e-10)
    uiwait (errordlg ('Axis limits are identical.'));
    status = 0;
    return
  end
  
  % NaNs will be filled in later after the values are computed...
  DEPLOT.x = [x_min,x_max,NaN,NaN];
  DEPLOT.y = [NaN,NaN,y_min,y_max];
return

function [X,Y,action] = get_point (marker,flag)
  global DEPLOT
  if (nargin < 1), marker = 'r+'; end
  if (nargin < 2), flag = 'data'; end
  h = [];
  action = 'continue';
  while(~strcmp(action,'point') & ~strcmp(action,'quit') & ~strcmp(action,'new'))
    [X,Y,button] = ginput(1);
    if (ishandle (h))
      delete (h);
    end
    switch (button)
     case double('q')			% quit digitizing
      unzoom;
      action = 'quit';
      X = NaN; Y = NaN;
     case double('z')			% zoom in
      limits = axis;
      zoom_factor = 7;
      dX = (limits(2)-limits(1))/zoom_factor;
      dY = (limits(4)-limits(3))/zoom_factor;
      % sanity
      dX = max(1,dX);	dY = max(1,dY);
      zoom_limits = [round(X-dX),round(X+dX),round(Y-dY),round(Y+dY)];
      if (any (zoom_limits ~= limits))
	DEPLOT.zoom_level = DEPLOT.zoom_level + 1;
	DEPLOT.zoom_limits{DEPLOT.zoom_level} = zoom_limits;
	axis (zoom_limits);
      end
      action = 'zoom';
      X = NaN; Y = NaN;
     case double ('u')			% unzoom one level
      if (DEPLOT.zoom_level > 1)
	DEPLOT.zoom_level = DEPLOT.zoom_level - 1;
	axis (DEPLOT.zoom_limits{DEPLOT.zoom_level});
      end
      action = 'unzoom';
      X = NaN; Y = NaN;
     case double ('U')			% unzoom all the way
      unzoom;
      action = 'unzoom';
      X = NaN; Y = NaN;
     case double ('n')			% new curve
      unzoom;
      action = 'new';
      DEPLOT.curve_id = DEPLOT.curve_id + 1;
      X = NaN; Y = NaN;
     case {1,2,3}			% select a point with mouse
      action = 'point';
      plot(X,Y,marker);			% mark selected point
     case double ('s')			% show coordinates
      action = 'show';
      h = show_coordinates(X,Y);
     otherwise				% shouldn't happen (ha, ha, ha)
      action = 'continue';
      X = NaN; Y = NaN;
    end
  end
  if (ishandle(h))
    delete (h);
  end
  % 2018.06.18 CB: commenting below to stop zoom-out after clicking a point
%   if (~strcmp(action,'point') | strcmp(flag,'axis'))
%     % unzoom all the way...
%     unzoom;
%   else
%     % unzoom one level...
%     if (DEPLOT.zoom_level > 1)
%       DEPLOT.zoom_level = DEPLOT.zoom_level - 1;
%       axis (DEPLOT.zoom_limits{DEPLOT.zoom_level});
%     end
%   end
return

function status = load_image
  global DEPLOT;
  status = 1;
  [file,path] = uigetfile('*.*', 'Open Scanned File:');
  if (path == 0)
    status = 0;
  else
    % Read in file and plot
    try
      [Plot,map]=imread([path,file]);
      figure (gcbf); hold off;
      DEPLOT.xx = image(Plot);
      unzoom;
      if (isempty(map))
	map = gray;			% a good? guess
      end
      colormap(map);
      hold on;
    catch
      uiwait (errordlg ('Error reading image file.'));
      status = 0;
    end
  end
return

function status = save_data
  global DEPLOT;
  status = 1;
  [file,path] = uiputfile('*.txt', 'Save Data As:');
  if (path == 0)
    status = 0;
    uiwait (errordlg ('Error specifying output file.'));
  else
    data = DEPLOT.data;
    % write out header and data to file...
    fp = fopen ([path,file],'w');
    if (fp < 0)
      status = 0;
      uiwait (errordlg ('Cannot open output file.'));
      return
    end
    fprintf (fp,'%% Deplotted %s\n',datestr(now,0));
    fprintf (fp,'%% x\ty\tcurve-id\n');
    [rows,cols] = size (data);
    for r = 1:rows
      % we do this by hand because 'save -ascii -append' doesn't work
      fprintf (fp,'%.6g\t%.6g\t%d\n',data(r,1),data(r,2),data(r,3));
    end
    fclose (fp);
    % save successful, so disable the button...
    set (findobj(gcbf,'Tag','SaveButton'),'Enable','off');
  end
return

function status = digitize_axes
  global DEPLOT
  status = 1;
  
  % CAS Jan 5 2006
  % For some reason, I now need to give the image the explicit mouse focus
  % by clicking on it (is there some way to do this programatically?)...
  title ('Click mouse anywhere on the image.','Color',[1,0,0]);
  [X,Y,button] = ginput(1);		% ignore

  if (DEPLOT.explicit_origin)
    title ('Click on origin.','Color',[1,0,0]);
    [X1,Y1,action] = get_point ('b+','axis');
    X3 = X1; Y3 = Y1;
  else
    title ('Click on X-axis minimum.','Color',[1,0,0]);
    [X1,Y1,action] = get_point ('b+','axis');
  end
  if (isnan(X1) | isnan(Y1))
    uiwait (errordlg ('Error digitizing point.'));
    status = 0;
    return
  end

  title ('Click on X-axis maximum.','Color',[1,0,0]);
  [X2,Y2,action] = get_point ('b+','axis');
  if (isnan(X2) | isnan(Y2))
    uiwait (errordlg ('Error digitizing point.'));
    status = 0;
    return
  end

  if (~DEPLOT.explicit_origin)
    title ('Click on Y-axis minimum.','Color',[1,0,0]);
    [X3,Y3,action] = get_point ('b+','axis');
    if (isnan(X3) | isnan(Y3))
      uiwait (errordlg ('Error digitizing point.'));
      status = 0;
      return
    end
  end
  title ('Click on Y-axis maximum.','Color',[1,0,0]);
  [X4,Y4,action] = get_point ('b+','axis');
  if (isnan(X4) | isnan(Y4))
    uiwait (errordlg ('Error digitizing point.'));
    status = 0;
    return
  end
  
  DEPLOT.X = [X1,X2,X3,X4];
  DEPLOT.Y = [Y1,Y2,Y3,Y4];
return

% compute transformation matrix...
function status = find_matrix
  global DEPLOT;
  status = 1;

  [x1,x2,x3,x4]=explode(DEPLOT.x);
  [y1,y2,y3,y4]=explode(DEPLOT.y);
  [X1,X2,X3,X4]=explode(DEPLOT.X);
  [Y1,Y2,Y3,Y4]=explode(DEPLOT.Y);

  dx = x2 - x1;
  dy = y4 - y3;

  a_11 = (X2 - X1)/dx;
  a_21 = (Y2 - Y1)/dx;
  a_22 = (Y4 - Y3)/dy;
  a_12 = (X4 - X3)/dy;

  A = [a_11, a_12; a_21, a_22];
  detA = det(A);

  if (abs(detA)<1e-10)
    uiwait (errordlg ('Points do not span a coordinate frame.'));
    status = 0;
    return;
  end

  B = inv(A)*[X1-X3;Y1-Y3];
  x3 = x1 - B(1); 
  x4 = x3; 				% by definition of axis
  y1 = y3 + B(2); 
  y2 = y3;				% by definition of axis

  Xt = [X1;Y1] - A*[x1;y1];

  if (DEPLOT.draw_axes | DEPLOT.force_perp_axes)
    % find origin...
    B = A*[x3;y1] + Xt;
    X0 = B(1); Y0 = B(2);
  end

  % if forcing perp axes, we'll draw them on the second time through
  if (DEPLOT.draw_axes & ~DEPLOT.force_perp_axes)
    plot (X0,Y0,'bo');
    if (norm([X1-X0,Y1-Y0])>norm([X2-X0,Y2-Y0]))
      X = X1; Y = Y1;
    else
      X = X2; Y = Y2;
    end
    plot ([X0,X],[Y0,Y],'b:');
    if (norm([X3-X0,Y3-Y0])>norm([X4-X0,Y4-Y0]))
      X = X3; Y = Y3;
    else
      X = X4; Y = Y4;
    end
    plot ([X0,X],[Y0,Y],'b:');
  end

  if (DEPLOT.force_perp_axes)
    DEPLOT.force_perp_axes = 0;		% do this only once
    v1 = [X2-X0;Y2-Y0];
    v2 = [X4-X0;Y4-Y0];
    
    % compute angle between the axes...
    theta = acos (dot(v1,v2)/(norm(v1)*norm(v2)));
    % compute angle by which to rotate each axis...
    phi = -(pi/2 - theta)/2;
    % in normal case rotate x axis by -phi
    if ((X2>X0 & Y4>Y0) | (X2<X0 & Y4<Y0))
      angle = -phi;
    else
      angle = phi;
    end
    
    R = [cos(angle), sin(angle); 
	 -sin(angle), cos(angle)];
    
    % rotate x axis
    pt = [X0;Y0] + R*[X1-X0;Y1-Y0];
    X1 = pt(1); Y1 = pt(2);

    pt = [X0;Y0] + R*[X2-X0;Y2-Y0];
    X2 = pt(1); Y2 = pt(2);  

    % rotate y axis in opposite direction
    angle = -angle;
    R = [cos(angle), sin(angle); 
	 -sin(angle), cos(angle)];
    
    pt = [X0;Y0] + R*[X3-X0;Y3-Y0];
    X3 = pt(1); Y3 = pt(2);  

    pt = [X0;Y0] + R*[X4-X0;Y4-Y0];
    X4 = pt(1); Y4 = pt(2);  

    DEPLOT.X = [X1,X2,X3,X4];
    DEPLOT.Y = [Y1,Y2,Y3,Y4];
    if (~find_matrix)	% recompute the transformation matrix
      status = 0;
      return;
    end
  end

  % save transformation matrix...
  DEPLOT.invA = inv(A);
  DEPLOT.Xt = Xt;
  DEPLOT.transform_defined = 1;
return

function data = digitize_points
  global DEPLOT
  point_styles = ['o','x','*','+','s','d','v','^','<','>','p','h'];
  DEPLOT.curve_id = 0;
  data = [];
  action = 'ok';
  lastpt = [];
  title ('Click on points (z=zoom,u=unzoom,n=new,q=quit).','Color',[1,0,0]);
  while (~strcmp(action,'quit'))
    % cycle through the available point styles...
    style = point_styles(1+rem(DEPLOT.curve_id,length(point_styles)));
    [X,Y,action] = get_point (['r',style],'data');
    if (strcmp(action,'point'))		% we got a point!
      [x,y] = transform (X,Y);
      data = [data; [x,y,DEPLOT.curve_id]];
      % draw connecting line if previous point belongs to the same curve
      if (~isempty(lastpt) & (lastpt(3) == DEPLOT.curve_id))
	plot ([lastpt(1),X],[lastpt(2),Y],'r-');
      end
      lastpt = [X,Y,DEPLOT.curve_id];
    end
  end
return

function [x,y] = transform (X,Y);
% transform the point X,Y in image coordinates to
% the graph coordinates...
  global DEPLOT

  invA = DEPLOT.invA;
  Xt = DEPLOT.Xt;

  X_pt = [X;Y];
  x_pt = invA*(X_pt-Xt);

  x = x_pt(1);
  y = x_pt(2);

  if (DEPLOT.x_log)
    x = exp (x);
  end
  if (DEPLOT.y_log)
    y = exp (y);
  end
return

% enable/disable all buttons and panels...
% (state is either 'on' or 'off')
function enable_buttons (state)
  global DEPLOT;
  set (findobj(gcbf,'Tag','XLogCheckbox'),'Enable',state);
  set (findobj(gcbf,'Tag','YLogCheckbox'),'Enable',state);
  set (findobj(gcbf,'Tag','OriginCheckbox'),'Enable',state);
  set (findobj(gcbf,'Tag','DrawAxesCheckbox'),'Enable',state);
  set (findobj(gcbf,'Tag','PerpAxesCheckbox'),'Enable',state);
  set (findobj(gcbf,'Tag','XminEditText'),'Enable',state);
  set (findobj(gcbf,'Tag','XmaxEditText'),'Enable',state);
  set (findobj(gcbf,'Tag','YminEditText'),'Enable',state);
  set (findobj(gcbf,'Tag','YmaxEditText'),'Enable',state);
  set (findobj(gcbf,'Tag','LoadButton'),'Enable',state);
  set (findobj(gcbf,'Tag','DigitizeButton'),'Enable',state);
  set (findobj(gcbf,'Tag','QuitButton'),'Enable',state);
  set (findobj(gcbf,'Tag','HelpButton'),'Enable',state);
  set (findobj(gcbf,'Tag','SaveButton'),'Enable',state);
  if (strcmp(state,'on') & isempty(DEPLOT.data))
    % disable button since there's nothing to save...
    set (findobj(gcbf,'Tag','SaveButton'),'Enable','off');
  end
return

function [varargout] = explode(X)
% function [x1,x2,...] = explode (X)
% Distributes (explodes) x column-wise into the output variables.
% C.A.Shera
  [m,n]=size(X);
  if (nargout ~= n)
    error ('explode(X): nargout != cols(X)');
  end

  for i = 1:nargout
    eval([ 'varargout(', int2str(i) ,')={X(:,', int2str(i), ')};' ]);
  end
return

function unzoom;
  global DEPLOT;
  axis image;
  DEPLOT.zoom_level = 1;
  clear DEPLOT.zoom_limits;
  DEPLOT.zoom_limits{1} = axis;
return

function h = show_coordinates (X,Y)
  global DEPLOT;
  if (DEPLOT.transform_defined)
    [x,y] = transform(X,Y);
  else
    x = X; y = Y;
  end
  h = text (X,Y,['(',num2str(x,'%.4g'),',',num2str(y,'%.4g'),')'], ...
	    'Color',[1,0,0],'FontSize',10);
return
