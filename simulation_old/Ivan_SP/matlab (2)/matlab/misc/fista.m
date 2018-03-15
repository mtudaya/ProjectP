function [x,J] = fista(y,H,Ht,lambda,alpha,Nit,x_init)
% [x, J] = fista(y,H,Ht,lambda,alpha,Nit)
%
% l1 regularization with 'fast iterated soft-thresholding algorithm'
% by Beck and Teboulle,
% Minimizes J(x) = ||y-H*x||_2^2 + lambda*||x||_1
%  INPUT
%   y - observed signal
%   H - function handle for operator
%   Ht - function handle for operator conjugate transpose
%   lambda - regularization parameter
%   alpha - need alpha > max(eig(H'*H))
%   Nit - Number of iterations
% OUTPUT
%   x - result of deconvolution
%   J - objective function
%
% To specificy initialization:
% [x,J] = fista(y,H,Ht,lambda,alpha,Nit,x_init)

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu

J = zeros(1, Nit);       % Objective function
if nargin < 7 
    x = 0*Ht(y);             % Initialize x
else
    x = x_init;
end
g = x;
T = lambda/(2*alpha);
t = zeros(1,Nit+1);
t(1) = 1;
for k = 1:Nit
    Hx = H(x);
    J(k) = sum(abs(Hx(:)-y(:)).^2) + lambda*sum(abs(x(:)));
    
    Hg = H(g);
    x_old = x;
    x = soft_fun(g + (Ht((y - Hg)))/alpha, T);
    
    t(k+1) = (1+sqrt(1+4*t(k)^2))/2;
    g = x + (t(k)-1)/t(k+1) * (x - x_old);            
    
end


% --- local functions ----

function y = soft_fun(x,T)
% Soft-threshold function
% This version works for real and complex data.
% y = soft_fun(x,T)
% x : input data
% T : threshold

g = max(abs(x)-T,0);
y = g./(g+T) .* x;

