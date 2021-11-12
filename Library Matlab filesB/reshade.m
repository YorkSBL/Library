function y= reshade(x,In);


if (In.cmin>= In.cmax), error('Inapprpriate bounds (cmin>=cmax)'); end

if (x>=In.cmax)    
    c=0;  
elseif (x<=In.cmin)    
    c=In.alpha;
else
    %c= 0.5*In.alpha*(1+(In.cmax+In.cmin)/(In.cmax-In.cmin)) - In.alpha*x/(In.cmax-In.cmin);
    c=  In.alpha/(In.cmax-In.cmin) * (In.cmax-x);

end


y=c;

end