function Fit = loessfit (xi,yi,Opt)
% function Fit = loessfit (xi,yi,Opt)
% Returns the loess fit y(x) to data points yi(xi) 
% with options specified in the struct Opt.
% The input arguments are:  
%   xi   : x values of data points
%   yi   : y values of data points
%   Opt  : struct of options with fields [defaults in brackets]
%    Opt.x [Nx pts spanning xi] : x values at which to compute the fit
%    Opt.Nx [100]               : number of default output points
%    Opt.span [0.2]             : size of moving window 
%    Opt.spanmeaning [standard] : how to interpret Opt.span [point fraction]
%    Opt.autospan [F]           : determine span using cross-validation
%    Opt.kfold [10]             : folds k for autospan cross-validation
%    Opt.lambda [1]             : order of the fitting polynomial
%    Opt.weights []             : weights to use in fitting (e.g., 1/sigma^2)
%    Opt.groupids []            : group ids of data points (see Opt.resample)
%    Opt.robust [0]             : number of times to iterate for robustness
%    Opt.collapse [F]           : replace points at same xi with median?
%    Opt.dither [0]             : percent by which to dither the xi
%    Opt.uselogx [F]            : log transform xi values before fitting
%    Opt.uselogy [F]            : log transform yi values before fitting
%    Opt.removeInfs [T]         : remove Infs from input?
%    Opt.removeNaNs [T]         : remove NaNs from input?
%    Opt.verbose [T]            : spew out possibly helpful information?
%    Opt.resample.data [F]      : resample data to estimate fit stddev?
%    Opt.resample.span [F]      : resample loess span? 
%    Opt.resample.spanrange     : range to vary span over [span*[2/3,3/2]]
%    Opt.resample.N [100]       : number of resamples to compute
%    Opt.resample.CI [95]       : % confidence intervals to compute
%    Opt.resample.bygroupid [F] : resample by selecting from group ids
%    Opt.resample.returnall [F] : return all N resampled trends
% The output structure Fit contains:
%   Fit.x                       : x values for fit 
%   Fit.y                       : the loess fit
%   Fit.dy                      : standard deviation of fit (if resampled)
%   Fit.CI.Upper                : upper bound of CI of fit (if resampled)
%   Fit.CI.Lower                : lower bound of CI of fit (if resampled)
%   Fit.CI.All                  : all trends from resampled data (if resampled)
%   Fit.DataUsed.xi             : xi actually used    
%   Fit.DataUsed.yi             : xi actually used
%   Fit.Opt                     : the options used for this fit
%
% The smoothing parameter Opt.span specifies the size of the moving
% window. The interpretation of span is given by the string spanmeaning.
% The options for Opt.spanmeaning are:
%   'standard'         : fraction of points, numel(xi) [aka 'point fraction']
%   'range fraction'   : fraction of xi range, max(xi)-min(xi)
%   'interval'         : length of xi interval (a linear distance)
%   'octaves'          : number of xi octaves (a log distance)
% The default and conventional meaning of span is obtained with 'standard';
% span then gives the window size as a fraction of the total 
% number of data points. Typically, span<<1 (so that the fit is *local*).
% Span can be either a constant (scalar) or an array (same size as
% Opt.x) that varies with x. Opt.span(i) is then the value of
% span to use at Opt.x(i)
%
% The fitting uses a polynomial fit of order Opt.lambda. 
% Note that the data (xi,yi) do NOT need to be pre-sorted by xi.
%
% The robust option iterates to compute a robust loess fit.
% First, an initial (non-robust) loess fit is performed. The
% weights are then modified (based on the residuals to reduce
% the influence of outliers) and a new fit is performed. 
% Ideally, the cycle of fitting and adjusting weights based on 
% residuals is iterated until the fit converges. 
% Here, the fit is iterated a total of robust times; 
% in practice robust=1 often suffices. If there are lots of
% data points, a robust fit can take awhile.
%
% The collapse option collapses all data points at the same
% value of xi and replaces the multiple yi by their median value.
% (We use the median not the mean in the hope that it'll be
% less sensitive to outliers.) Collapsing helps maintain 
% the 'local' character of the fit in cases where many xi have 
% multiple corresponding yi.
%  
% The value of dither specifies the percent by which to dither
% the xi randomly using a Gaussian distribution. Dithering 
% is an alternative and generally better solution to the 
% 'duplicate xi problem' addressed by collapsing. If dither is 
% set, collapse is unset.  
%
% The resample options allow one to estimate the uncertainty
% in the fit by resampling either the data (xi,yi) with replacemet
% and/or the value of span over a specified range and recomputing
% the fit N times. The output variable dy gives the resulting
% standard deviation of the fit at each x value.
%  
% Reference: 
%  Visualizing Data, William S. Cleveland
%  Hobart Press, 1993, pg 100ff
%
%  C.A. Shera
%
  
  if (nargin<3), Opt = []; end
  
  if (~isusedfield(Opt,'x')),              Opt.x = []; end
  if (~isusedfield(Opt,'Nx')),             Opt.Nx = 100; end
  if (~isusedfield(Opt,'span')),           Opt.span = 0.2; end
  if (~isusedfield(Opt,'spanmeaning')),    Opt.spanmeaning = 'standard'; end  
  if (~isusedfield(Opt,'autospan')),       Opt.autospan = false; end  
  if (~isusedfield(Opt,'kfold')),          Opt.kfold = 10; end  
  if (~isusedfield(Opt,'lambda')),         Opt.lambda = 1; end
  if (~isusedfield(Opt,'weights')),        Opt.weights = []; end
  if (~isusedfield(Opt,'groupids')),       Opt.groupids = []; end
  if (~isusedfield(Opt,'robust')),         Opt.robust = 0; end
  if (~isusedfield(Opt,'collapse')),       Opt.collapse = false; end
  if (~isusedfield(Opt,'dither')),         Opt.dither = 0; end
  if (~isusedfield(Opt,'uselogx')),        Opt.uselogx = false; end
  if (~isusedfield(Opt,'uselogy')),        Opt.uselogy = false; end
  if (~isusedfield(Opt,'removeInfs')),     Opt.removeInfs = true; end
  if (~isusedfield(Opt,'removeNaNs')),     Opt.removeNaNs = true; end
  if (~isusedfield(Opt,'verbose')),        Opt.verbose = true; end
  if (~isusedfield(Opt,'resample.data')),  Opt.resample.data = false; end
  if (~isusedfield(Opt,'resample.span')),  Opt.resample.span = false; end
  if (~isusedfield(Opt,'resample.N')),     Opt.resample.N = 100; end
  if (~isusedfield(Opt,'resample.CI')),    Opt.resample.CI = 95; end
  if (~isusedfield(Opt,'resample.bygroupid')), 
    Opt.resample.bygroupid = false; 
  end
  if (~isusedfield(Opt,'resample.spanrange'))
    Opt.resample.spanrange = Opt.span(:)*[2/3,3/2]; 
  end
  if (~isusedfield(Opt,'resample.returnall')), 
    Opt.resample.returnall = false; 
  end
  
  % for backward compatibility...
  if (isusedfield(Opt,'alpha'))
    Opt.span = Opt.alpha;
  end
  if (isusedfield(Opt,'alphameaning')), 
    Opt.spanmeaning = Opt.alphameaning;
  end
  if (isusedfield(Opt,'resample.alpha'))
    Opt.resample.span = Opt.resample.alpha;
  end
  if (isusedfield(Opt,'resample.alpharange'))
    Opt.resample.spanrange = Opt.resample.alpharange;
  end

  if (strcmp(Opt.spanmeaning,'point fraction'))
    Opt.spanmeaning = 'standard';      % synonyms
  end
  
  resampling = (Opt.resample.data || Opt.resample.span);
  
  if (Opt.resample.bygroupid && ~resampling)
    Opt.resample.bygroupid = false;
  end
  
  span_is_array = false;
  if (numel(Opt.span) > 1)
    span_is_array = true;
    if (numel(Opt.span) ~= numel(Opt.x))
      error ('size of span array does not match Opt.x');
    end
    if (Opt.resample.span)
      error ('resampling of vector span not supported');
    end
  end
  
  if (Opt.autospan)
    if (span_is_array)
      error ('autospan requires scalar span');
    end
    if (~strcmp(Opt.spanmeaning,'standard'))
      error ('autospan requires standard spanmeaning');
    end
    if (Opt.resample.span)
      error ('cannot resample autospan');
    end
  end
  
  if (Opt.dither), Opt.collapse = false; end

  W = abs(Opt.weights);                 % must be non-negative
  if (isempty(W))
    W = ones(size(yi));
  end

  % do sanity control on the robust flag since
  % we later use it to control the iteration... 
  Opt.robust = round (abs (Opt.robust));
  
  % do some minimal sanity checking...
  if (Opt.resample.span & ~inrange(Opt.span,Opt.resample.spanrange))
    error ('loess span not in resampling range');
  end
  if (resampling && Opt.resample.N<2)
    Opt.resample.data = false; Opt.resample.span = false;
    if (Opt.verbose)
      disp('loessfit: N too small to resample');
    end
  end
  if (Opt.resample.bygroupid)
    gids = Opt.groupids(:);  % make columnar
    if (numel(gids) ~= numel(xi))
      Opt.resample.bygroupid = false;
      if (Opt.verbose)
        disp('loessfit: group ids do not match data');
      end
    end
  end
  if (any(isinf(xi) | isinf(yi)))
    if (Opt.removeInfs)
      disp ('Removing Infs in input data');
      ok = ~(isinf(xi) | isinf(yi));
      xi = xi(ok); yi = yi(ok); W = W(ok); 
      if (Opt.resample.bygroupid), gids = gids(ok); end
    else
      error ('Infs in input data');
    end
  end
  if (any(isnan(xi) | isnan(yi)))
    if (Opt.removeNaNs)
      disp ('Removing NaNs in input data');
      ok = ~(isnan(xi) | isnan(yi));
      xi = xi(ok); yi = yi(ok); W = W(ok);
      if (Opt.resample.bygroupid), gids = gids(ok); end
    else
      error ('NaNs in input data');
    end
  end
  
  % collapse repeated values...
  if (Opt.collapse)
    [sxi,I,J] = unique (xi);
    syi = zeros(size(sxi));
    wgt = zeros(size(sxi));
    for k = 1:numel(sxi)
      % FIXME: Note that the computation of the mean or median 
      % should take the weights into account and new weights 
      % should be computed using the variance.
      avg_me = find (J==k);
      syi(k) = median (yi(avg_me));     % should use weights
      wgt(k) = mean (W(avg_me));        % coarse
    end
    xi = sxi;
    yi = syi;
    W = wgt;
  end

  % possibly do log transforms...
  xio = xi;                             % save original
  if (Opt.uselogx)
    ok = (xi>0);
    if (~all(ok))
      if (Opt.verbose), disp('loessfit: excluding xi<0 for log'); end
      xi = xi(ok); yi = yi(ok); W = W(ok); xio = xio(ok);
      if (Opt.resample.bygroupid), gids = gids(ok); end
    end
    xi = log(xi);
    
    if (~isempty(Opt.x))
      ok = (Opt.x>0);
      Opt.x = log(Opt.x(ok));
    end
  end
  if (Opt.uselogy)
    ok = (yi>0);
    if (~all(ok))
      if (Opt.verbose), disp('loessfit: excluding yi<0 for log'); end
      xi = xi(ok); yi = yi(ok); W = W(ok);
      if (Opt.resample.bygroupid), gids = gids(ok); end
    end
    yi = log(yi);
  end
  
  % choose x values at which to compute the fit...
  if (isempty(Opt.x))
    % xi is already in log domain if Opt.uselogx
    Opt.x = linscale (min(xi),max(xi),Opt.Nx,isrow(xi));
  end
  if (Opt.robust > 0)
    xR = xi;
    interpolate_on_robust = false;
    if (numel(xi) > Opt.Nx)             % use fewer points then interpolate
      % xi is already in log domain if Opt.uselogx
      xR = linscale (min(xi),max(xi),Opt.Nx,isrow(xi));
      interpolate_on_robust = true;
    end
  end

  % possibly dither the xi...
  if (Opt.dither ~= 0)
    xi_orig = xi;                       % save originals
    runi = -1 + 2*rand(size(xi));       % [-1,1]
    xi = xi_orig .* (1 + Opt.dither*runi/100);
  end
  
  % possibly reinterpret span, converting to point fraction when necessary...
  if (span_is_array || strcmp(Opt.spanmeaning,'standard'))
    % do nothing (array-valued spans are assumed to be point fractions)
  else
    Nxi = numel(xi);
    spans = zeros(size(Opt.x));
    xrng = Opt.span;
    minxio = min(xio);
    maxxio = max(xio);
    if (strcmp(Opt.spanmeaning,'range fraction'))
      xrng = xrng*abs(maxxio-minxio);
    elseif (strcmp(Opt.spanmeaning,'octaves'))
      log2minxio = log2(minxio); 
      log2maxxio = log2(maxxio);
      log2xio = log2(xio);
    end
    for i = 1:numel(Opt.x)
      thex = Opt.x(i);
      if (Opt.uselogx) 
        thex = exp(thex);               % undo log transform
      end
      switch (Opt.spanmeaning)
        case {'interval','range fraction'}
          x1 = min(max(minxio,thex-xrng/2),maxxio-xrng);
          % compute point fractions...
          spans(i) = sum(inrange(xio,[x1,x1+xrng]))/Nxi;
        case {'octaves'}
          x1 = min(max(log2minxio,log2(thex)-xrng/2),log2maxxio-xrng);
          % compute point fractions...
          spans(i) = sum(inrange(log2xio,[x1,x1+xrng]))/Nxi;
      end
    end
    % save original settings for the record...
    Opt.span0 = Opt.span;
    Opt.span0meaning = Opt.spanmeaning;
    % new settings...
    Opt.span = spans;
    Opt.spanmeaning = 'standard';
    span_is_array = true;
  end

  % possibly use autospan...
  if (Opt.autospan)
    OptR = Opt; 
    % modify and/or turn off options we don't want...
    OptR.robust = 0; 
    OptR.dither = 0;                    % points already dithered
    OptR.autospan = false;
    OptR.spanmeaning = 'standard';
    OptR.uselogx = false;               % data already transformed
    OptR.uselogy = false;               % data already transformed
    OptR.resample.span = false;
    OptR.resample.data = false;
    OptR.resample.bygroupid = false;
    OprR.weights = W;;

    % See http://blogs.mathworks.com/loren/2011/01/13/data-driven-fitting/
    cp = cvpartition(numel(xi),'Kfold',Opt.kfold);
    Opt.span = fminbnd(@(x) myloessfitfun(x,[xi(:),yi(:)],OptR,cp), ...
                       0.01,0.99,optimset('TolX',0.01));
    % now that we know span, set autospan to false...
    Opt.autospan = false;
    Opt.autospan0 = true;
    if (Opt.verbose)
      disp(['loessfit: autospan gives span of ',num2str(Opt.span)]);
      if (isnan(Opt.span))
        error('loessfit: autospan returned NaN');
      end
    end
  end
 
  % possibly iterate to obtain robust fits at the points xi...
  if (Opt.robust > 0)
    original_weights = W;                 % save for later
    OptR = Opt; 
    % for speed, we compute fit at xR, then interpolate to xi...
    OptR.x = xR; 
    if (span_is_array)
      % must interpolate to xR
      OptR.span = interp1(Opt.x,Opt.span,xR);
    end
    
    % modify and/or turn off options we don't want...
    OptR.robust = 0; 
    OptR.dither = 0;                    % points already dithered
    OptR.spanmeaning = 'standard';
    OptR.uselogx = false;               % data already transformed
    OptR.uselogy = false;               % data already transformed
    OptR.resample.span = false;
    OptR.resample.data = false;
    OptR.resample.bygroupid = false;

    robust = Opt.robust;
    while (robust > 0)
      robust = robust - 1;
    
      % do a normal loess fit at points xi (not x) using the current weights
      OprR.weights = W;;
      FitR = loessfit(xi,yi,OptR);

      % possibly interpolate...
      yr = FitR.y;
      if (interpolate_on_robust)
        % extrap is necessary because of dithering...
        yr = interp1(FitR.x,FitR.y,xi,'linear','extrap');
      end
      
      % compute residuals from fit...
      res = yi - yr;
    
      % Compute bisquare robustness weighting function using the residuals...
      % Outliers (points with large residuals) receive a weighting near zero.
      mar = median (abs (res));		% median absolute residual
      u = res / (6*mar);

      % before iterating the fit, modify the original weights 
      % using the bisquare robustness weighting...
      W = original_weights .* bisquare (u);
    end
  end

  % allocate space for the fit...
  y = zeros (size(Opt.x));
  
  span = Opt.span;
  
  % after this loop, robust is always zero...
  % Loop over x (not xi) ...
  for i = 1:numel(Opt.x)
    % Compute the tricube weighting functions...
    % Points are weighted by their distance from x(i) using a 
    % variable window defined so that the qth most distant point 
    % has w=0. Except near the ends, the window will typically 
    % be roughly symmetric about x(i). 
    
    % NB: We no longer do this...
    % Note that when determining 
    % the qth most distant point, we count points at the same xi 
    % as one. This differs from Cleveland's description, but is
    % more robust (and sensible?) in certain pathological cases.
    
    if (span_is_array)
      span = Opt.span(i);
    end
      
    Delta = abs(xi-Opt.x(i));
    %    Delta_q = unique (Delta);       % sort array, removing duplicates
    Delta_q = sort (Delta);       % sort array
    n = numel(Delta_q);
    q = max(1,min (n, round (abs (span)*n)));
    
    % multiply in the tricube weighting...
    w = W .* tricube (Delta/Delta_q(q));
    
    % Locate the points to fit...
    fit_me = find (w>0);

    % take care of various pathological cases...
    switch (numel(fit_me))
     case {0}
      if (Opt.verbose), disp('loessfit: Empty window! Using NaN'); end
      y(i) = NaN;
      continue;
     
     case {1}
      if (Opt.verbose), disp('loessfit: One point in window! Skipping fit'); end
      y(i) = yi(fit_me);
      continue;
     
     otherwise
      lam = min(Opt.lambda,numel(fit_me)-1);
      if (lam ~= Opt.lambda && Opt.verbose)
	disp('loessfit: Too few points in window! Reducing lambda');
      end

      % Do the fit...
      if (lam == 1)			
	% MUCH faster than polyfitw in linear case
	[a,b] = linear_fit (xi(fit_me),yi(fit_me),1./sqrt(w(fit_me)));
	y(i) = a + b*Opt.x(i);
      else
	p = polyfitw (xi(fit_me),yi(fit_me),lam,w(fit_me));
	y(i) = polyval (p,Opt.x(i));
      end
    end
  end

  if (Opt.dither ~= 0)
    xi = xi_orig;                       % restore undithered values
  end

  % invert the log transforms...
  if (Opt.uselogx)
    xi = exp(xi); Opt.x = exp(Opt.x);
  end
  if (Opt.uselogy)
    yi = exp(yi); y = exp(y);
  end

  if (Opt.verbose && any(isinf(y)))
    disp('loessfit: Infs in fit');
  end
  if (Opt.verbose && any(isnan(y)))
    disp('loessfit: NaNs in fit');
  end

  % store it all away...
  Fit.x = Opt.x;
  Fit.y = y;
  Fit.dy = [];                          % default value
  Fit.CI.Upper = [];                    % default value
  Fit.CI.Lower = [];                    % default value
  Fit.DataUsed.xi = xi;
  Fit.DataUsed.yi = yi;
  
  % possibly bootstrap some fits to compute std dev...
  if (resampling)
    yBoot = zeros(numel(Opt.x),Opt.resample.N);
    if (Opt.resample.bygroupid) 
      nGrps = numel(unique(gids));
    end
    if (Opt.verbose), fprintf('loessfit: Resampling 0'); end
    for n=1:Opt.resample.N
      if (Opt.verbose) 
        fprintf(qcolon(mod(n,5)==0,{num2str(n),'.'}));
        if (n~=Opt.resample.N & mod(n,50)==0), fprintf('\n'); end
      end
      if (Opt.resample.data)
        % random sample w/ replacement
        if (Opt.resample.bygroupid)
          % sample from groups
          J = floor(nGrps*rand(1,nGrps)+1);
          I = [];
          for id=J
            I = [I;find(gids==id)];   % gids is column vector
          end
        else
          % sample from all data
          I = floor(numel(xi)*rand(size(xi))+1);  
        end
      else
        I = 1:numel(xi);
      end
      if (Opt.resample.span)
        % random span in range
        span = Opt.resample.spanrange(1) + rand*diff(Opt.resample.spanrange);
      else
        span = Opt.span;
      end
        
      % modify and/or turn off options we don't want...
      OptB = Opt;
      OptB.span = span;
      OptB.resample.data = false;
      OptB.resample.span = false;
      OptB.resample.bygroupid = false;
      FitB = loessfit(xi(I),yi(I),OptB);
      yBoot(:,n) = FitB.y(:);
    end
    if (Opt.verbose), fprintf('\n'); end

    % do sanity checking...
    if (any(isnan(yBoot) | isinf(yBoot)))
      if (Opt.verbose) 
        disp('loessfit: Inf or NaNs in resampled fits'); 
      end
      ok = find(isinf(yBoot))
      yBoot(ok) = NaN;
    end
    Fit.dy = nanstd(yBoot');
    dp = (100-Opt.resample.CI)/2;
    %%Fit.CI.Upper = prctile(yBoot',100-dp);
    %%Fit.CI.Lower = prctile(yBoot',dp);
    if (iscolumn(Opt.x))
      Fit.dy = Fit.dy(:);
      %%Fit.CI.Upper = Fit.CI.Upper(:);
      %%Fit.CI.Lower = Fit.CI.Lower(:);
    end
    if (Opt.resample.returnall)
      Fit.CI.All = yBoot';
    end
  end
  
  Fit.Opt = Opt;
  return

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function T = tricube (u)
  u = abs (u);
  T = zeros (size (u));
  i = find (u<1);
  T(i) = (1-u(i).^3).^3;
  return

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function T = bisquare (u)
  u = abs (u);
  T = zeros (size (u));
  i = find (u<1);
  T(i) = (1-u(i).^2).^2;
  return



