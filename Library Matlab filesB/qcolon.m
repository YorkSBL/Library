function z = qcolon (bool, values)
% Implements a functional version of the ?: operator
% function z = qcolon (bool, values)
%    bool is a cell array (or array) of logical values
%    values is the corresponding cell array of values
% Examples:
%   qcolon({0}  ,{'dog','cat'})        returns 'cat'
%   qcolon({1,0},{'dog','cat','bear'}) returns 'dog'
%   qcolon({0,1},{'dog','cat','bear'}) returns 'cat'
%   qcolon({0,0},{'dog','cat','bear'}) returns 'bear'
% See also question_colon()
  
  n = numel(bool);
  if (numel(values) ~= n+1)
    error('argument mismatch');
  end

  if (~iscell(bool))
    bool = num2cell(bool);
  end
  if (~iscell(values))
    values = num2cell(values);
  end

  for k=1:n
    bool{k} = logical(bool{k});
  end
  bool{n+1} = logical(1);
  
  ok = find(logical(cell2mat(bool)));
  z = values{ok(1)};
  
