function [output] = maxF3(input)
%MAXF3 Summary of this function goes here
x = input(1);
y = input(2);

    term1 =  20.*exp((-0.2).*(x.^2+y.^2));
    term2 =  exp(cos(2.*pi.*x)+cos(2.*pi.*y));

    z = term1 + term2 ;
    
    output = z; 
end

