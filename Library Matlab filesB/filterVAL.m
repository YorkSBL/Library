function y = filterNAN(inputA,val)
% get rid of values that are equal to specified val    [C. Bergevin 8.23.11]
checkN= inputA == val;
mm=1; 
for nn=1:size(checkN,1)
  if checkN(nn)==0
    y(mm,:)=inputA(nn,:);
    mm=mm+1;
  end
end

return
