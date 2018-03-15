function [x, cost] = bp_salsa(y, A, At, p, mu, Nit, lambda)

% x = bp_salsa(y, A, At, p, mu, Nit)
%
% BASIS PURSUIT
% Minimize || x ||_1 such that y = A x
% where A At = p I
%
% INPUT
%   y     : data
%   A, At : function handles for A and its conj transpose
%   p     : Parseval constant
%   mu    : ADMM parameter (inverse step size)
%   Nit   : number of iterations
%
% OUTPUT
%   x     : solution to BP problem
%
% [x, cost] = bp_salsa(...) returns cost function history
%
% salsa_bp(..., lambda) minimizes ||lambda .* x||_1

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu
% March 2012

% The algorithm is a variant of SALSA (Afonso, Bioucas-Dias, Figueiredo,
% IEEE Trans Image Proc, 2010, p. 2345)

if nargin < 7
    lambda = 1;
end

% Initialization
x = At(y);
d = zeros(size(x));
cost = zeros(1,Nit);

for i = 1:Nit
    u = soft(x + d, 0.5*lambda/mu) - d;
    d = (1/p) * At(y - A(u));
    x = d + u;
    cost(i) = sum(lambda(:) .* abs(x(:)));
end
