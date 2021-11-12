function y = filterNAN(inputA)
% get rid of values that are NaN's    [C. Bergevin 7.16.08]
checkN=isnan(inputA);
mm=1; 
for nn=1:size(checkN,1)
  if checkN(nn)==0
    y(mm,:)=inputA(nn,:);
    mm=mm+1;
  end
end

return
