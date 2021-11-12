function z= myfun(x)
%z= AA./(BB-exp(-CC*x.^2)) + (rand(N,1)'*NA*AA);
z= x(1)./(x(2)-exp(-x(3)*x(6:end).^2)) + (rand(x(4),1)'*x(5)*x(1));