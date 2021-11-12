function z = eitheror (bool, true_value, false_value)
% Implements a functional version of the ?: operator
% Matlab is out of joint: O cursed spite, 
% That ever I was born to set it right!
% function z = eitheror (bool, true_value, false_value)

if (logical(bool))
  z = true_value;
else
  z = false_value;
end
