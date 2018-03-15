% AoAs/AoDs distribution
%AoAs ? and ? in the azimuth and elevation plane with reference to the respective specular direction
% modelled as  a zero-mean second order Gaussian Mixture Model (GMM)
%GMM(x) = a1* exp (-0.5*((x-x1_mu)/sigma1^2)))/(sqrt(2*pi))* sigma1 + a2* exp (-0.5*((x-x2_mu)/sigma1^2)))/(sqrt(2*pi))* sigma2
function[]=Angle()
a1=0.583;
sigma1=1.895;
x1_mu=0;
a2=0.417;
sigma2=5.487;
x2_mu=0;


x = -40:2:40; 
GMM= zeros(length(x),1);
for k=1:length(x)
xk = x(k);
GMM(k) = (a1* exp(-0.5*((xk-x1_mu)/sigma1)^2))/(sqrt(2*pi))* sigma1 + (a2* exp(-0.5*((xk-x2_mu)/sigma2)^2))/(sqrt(2*pi))* sigma2;
end
GMM=GMM./sum(GMM);
plot (x,GMM)

cdf = cumsum(GMM);
plot (x,cdf)

[cdf, mask] = unique(cdf);
x = x(mask);
% create an array of 2500 random numbers
randomValues = rand(1,400);

% inverse interpolation to achieve P(x) -> x projection of the random values
projection = interp1(cdf, x, randomValues);

hist(projection, 50);
%stem (projection);

%mean(projection)
