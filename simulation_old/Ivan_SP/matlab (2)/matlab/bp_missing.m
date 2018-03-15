function [c, cost] = bp_missing(y, F, Ft, p, s, mu, Nit, lambda)

% c = bp_missing(y, F, Ft, p, s, mu, Nit)
%
% MISSING data estimation using BP
% Minimize || c ||_1 such that y(s) = (F c)(s)
% where F Ft = p I
% The index vector s indicates known data.
%
% INPUT
%   y     : data
%   F, Ft : function handles for F and its conj transpose
%   p     : Parseval constant
%   s     : index vector for known data
%   mu    : ADMM parameter
%   Nit   : Number of iterations
%
% OUTPUT
%   c     : coefficients of estimated data
%
% [c, cost] = bp_missing(...) returns cost function history
%
% bp_missing(..., lambda) minimizes || lambda .* c ||_1

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu

% The algorithm is a variant of SALSA (Afonso, Bioucas-Dias, Figueiredo,
% IEEE Trans Image Proc, 2010, p. 2345)

if nargin < 8
    lambda = 1;
end

% Initialization
y(~s) = 0;
c = Ft(y)/p;
d = zeros(size(c));
cost = zeros(1,Nit);

for i = 1:Nit
    u = soft(c + d, 0.5*lambda/mu) - d;
    d = (1/p) * Ft(y - s.*F(u));
    c = d + u;
    cost(i) = sum(lambda(:) .* abs(c(:)));
end
