function isit = isusedfield(S,str)
% function isit = isusedfield(S,'field')
% Returns true if the specified field exists in the structure S and is non-empty

  isit = false;
  if (issubfield(S,str) & ~isempty(eval(strcat('S.',str,';'))))
    isit = true;
  end
  return
