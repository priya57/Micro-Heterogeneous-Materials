function S = Eshelby(a1,a2,mu)
S = zeros(6,6);

S(1,1) = 1/(2*(1-mu))*(((a2^2+2*a1*a2)/(a1+a2)^2)+((1-2*mu)*a2/(a1+a2)));
S(2,2) = 1/(2*(1-mu))*(((a2^2+2*a1*a2)/(a1+a2)^2)+((1-2*mu)*a1/(a1+a2)));
S(1,2) = 1/(2*(1-mu))*((a2^2/(a1+a2)^2)-((1-2*mu)*a2/(a1+a2)));
S(2,1) = 1/(2*(1-mu))*((a1^2/(a1+a2)^2)-((1-2*mu)*a1/(a1+a2)));
S(4,4) = 1/(2*(1-mu))*((a1^2+a2^2)/(2*(a1+a2)^2)+(1-2*mu)/2);
S(1,3) = mu/(2*(1-mu))*(2*a2/(a1+a2));
S(2,3) = mu/(2*(1-mu))*(2*a1/(a1+a2));
S(6,6) = a2/(2*(a1+a2));
S(5,5) = a1/(2*(a1+a2));
end
