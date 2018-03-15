% AoAs/AoDs distribution
%AoAs ? and ? in the azimuth and elevation plane with reference to the respective specular direction
% modelled as  a zero-mean second order Gaussian Mixture Model (GMM)
%GMM(x) = a1* exp (-0.5*((x-x1_mu)/sigma1^2)))/(sqrt(2*pi))* sigma1 + a2* exp (-0.5*((x-x2_mu)/sigma1^2)))/(sqrt(2*pi))* sigma2
function[projection]=Angle(a1,sigma1,x1_mu,a2,sigma2,x2_mu,data_points,angle_start,angle_end,angle_resolution)



x = angle_start:angle_resolution:angle_end; 
GMM= zeros(length(x),1);
for k=1:length(x)
xk = x(k);
GMM(k) = (a1* exp(-0.5*((xk-x1_mu)/sigma1)^2))/(sqrt(2*pi))* sigma1 + (a2* exp(-0.5*((xk-x2_mu)/sigma2)^2))/(sqrt(2*pi))* sigma2;
end
GMM=GMM./sum(GMM);
figure
plot (x,GMM)

cdf = cumsum(GMM);
%plot (x,cdf)

[cdf, mask] = unique(cdf);
x = x(mask);
% create an array of 2500 random numbers
randomValues = rand(1,data_points);

% inverse interpolation to achieve P(x) -> x projection of the random values
projection = interp1(cdf, x, randomValues);
projection(isnan(projection)) = []
%hist(projection, 50);
%stem (projection);

%mean(projection)

return
