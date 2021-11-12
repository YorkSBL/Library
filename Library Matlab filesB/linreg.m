function[a,b,m]=linreg(x,y,V)

clf

% data set
%x = [1 2 3 4 6 8 10 12 15 20 30 40 60];
%y = [118 58 38 31 22 13 14 12 6 7 6 2 1];


%V=1;  % error, assumed to be const. for all points (i.e. sigma_i)


% ************************************************************

%x=x.^-1
%y=y.^-1

N=size(x,2);

delta = N*sum(x.^2)-(sum(x))^2;

a=0;
%a = 1/delta * (sum(x.^2)*sum(y) - sum(x)*sum(x.*y));
b = 1/delta * (N*sum(x.*y) - sum(x)*sum(y));

sigmaA = 1/delta * (V^2*sum(x.^2));
sigmaB = 1/delta * (N*V^2);

chiS = sum( (1/V * (y-a-b*x)).^2);

S = N*(1/V)^2;
Sx = sum(x./V^2);
Sxx = sum((x.^2)/V^2);


rab = -Sx/(sqrt(S*Sxx));
Q = 1 - gammainc((N-2)/2,chiS/2);


plot(x,y,'*')
hold
plot([min(x):0.1:max(x)],b*[min(x):0.1:max(x)]+a)


temp=num2str(a);
temp2=num2str(b);
temp3=num2str(sigmaA);
temp4=num2str(sigmaB);
temp5=num2str(Q);
temp6=num2str(rab);
temp7=num2str(chiS);

text(max(x)-20,min(y)+0.4,['intercept is ',temp])
text(max(x)-20,min(y)+0.34,['slope is ',temp2])
text(max(x)-20,min(y)+0.28,['sigmaA is ',temp3])
text(max(x)-20,min(y)+0.22,['sigmaB is ',temp4])
text(max(x)-20,min(y)+0.16,['Q is ',temp5])
text(max(x)-20,min(y)+0.1,['Rab is ',temp6])
text(max(x)-20,min(y)+0.04,['chi^2 is ',temp7])




