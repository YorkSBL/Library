function tf = issubfield (S,fieldpath)
%ISSUBFIELD True if field is nested in structure array.
% ISSUBFIELD(S,'fieldpath') returns true if 'fieldpath' 
% is the name of a field or the path to a nested subfield
% in the structure array S. ISSUBFIELD extends the 
% functionality of isfield() to nested subfields 
%
% example: issubfield(A,'B.C') returns true if A.B.C exists
%
% Author: C.A.Shera
  
  tf = logical(0);
  if (isempty(fieldpath))
    return
  end

  [fields,nfields] = strsplit(fieldpath,'.');
  
  s1 = 'isfield(S';
  for n=1:nfields
    s2 = [',''',fields{n},''');'];
    ok = eval(strcat(s1,s2));
    if (~ok)
      return
    end
    s1 = strcat(s1,'.',fields{n});
  end

  tf = logical(1);
  return
  
  return