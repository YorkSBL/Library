% script to plot matrix as a vectorized "raster"

linewidth = 3;               % adjustable plotting parameter

% storeI is obtained by running Turbine simulation
if (~exist('storeI','var'))
  % make a random test case...
  storeI = rand(100);
  storeI(storeI>0.5) = 1;
  storeI(storeI<=0.5) = 0;
end

A = storeI;                  % a convenient synonym

% invert color scheme (switch 0s and 1s) to make it easier...
A = abs(A-1);
[nrows,ncols] = size(A);

x = 0:(ncols-1);
figure;
for m=1:nrows
  y = A(m,:);
  dy = [0,diff(y)];

  ups = find(dy>0);
  downs = find(dy<0);

  % deal with the ends...
  if (ups(1)>downs(1))
    ups = [1,ups];
  end
  if (ups(end)>downs(end))
    downs = [downs,ncols];
  end
  if (numel(ups) ~= numel(downs))
    warning('rethink this!');
  end
  
  for n=1:numel(ups)
    x1 = x(ups(n));
    x2 = x(downs(n));
    plot([x1,x2],[m,m],'k-','LineWidth',linewidth); hold on
  end
end

xlim([x(1),x(end)]);
ylim([1,nrows]);
