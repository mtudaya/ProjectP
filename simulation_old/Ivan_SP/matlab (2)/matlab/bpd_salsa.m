function [x, cost] = bpd_salsa(y, A, At, p, lambda, mu, Nit)

% x = bpd_salsa(y, A, At, p, lambda, mu, Nit)
%
% BASIS PURSUIT DENOISING
% Minimize ||y - A x||_2^2 + lambda * || x ||_1
% where A At = p I
%
% INPUT
%   y      : data
%   A, At  : function handles for A and its conj transpose
%   p      : Parseval constant
%   lambda : regularization parameter
%   mu     : ADMM parameter
%   Nit    : Number of iterations
%
% OUTPUT
%   x      : solution to BPD problem
%
% [x, cost] = bpd_salsa(...) returns cost function history

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu
% March 2012

% The program implements SALSA (Afonso, Bioucas-Dias, Figueiredo,
% IEEE Trans Image Proc, 2010, p. 2345)

if nargout > 1
    ComputeCost = true;
    cost = zeros(1,Nit);
else
    ComputeCost = false;
end    

x = At(y);
d = zeros(size(x));

for i = 1:Nit
    u = soft(x + d, 0.5*lambda/mu) - d;
    d = 1/(mu + p) * At(y - A(u));
    x = d + u;
    
    if ComputeCost
        residual = y - A(x);
        cost(i) = sum(abs(residual(:)).^2) + sum(abs(lambda * x(:))); 
    end
end

