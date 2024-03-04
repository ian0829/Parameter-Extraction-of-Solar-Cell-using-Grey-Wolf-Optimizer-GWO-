function [rmse, Ipv] =fobj(P)

global Vpv0 Ipv0

Rs      = P(1);
Rsh     = P(2);
Iph     = P(3);
Is      = P(4);
n       = P(5);
Vt      = 0.0257*36;

theta = ( ((Rsh * Rs)/(Rsh + Rs)) * (Is/1e9) * exp( (Rsh * Rs * (Iph + (Is/1e9)) + Rsh .* Vpv0)/(n * Vt * (Rsh + Rs))) ) / (n * Vt);
Ipv = (Rsh * (Iph + (Is/1e9)) - Vpv0)/(Rsh + Rs) - ( (n * Vt)/(Rs) ) * lambertw(theta);
Ipv(Ipv<0) = 0; 

 % write RMSE equation here 
rmse = sqrt(mean((Ipv0-Ipv).^2));
 
end
 