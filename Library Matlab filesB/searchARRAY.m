function array_position=searchARRAY(A,k)
% search array A for closest value k and return index  (CB 07.08)
[min_difference, array_position] = min(abs(A - k));
