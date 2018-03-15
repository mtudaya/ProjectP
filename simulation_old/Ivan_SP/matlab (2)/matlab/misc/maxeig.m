function p = maxeig(Afn,s,Nit)
% p = maxeig(Afn,s,Nit)
% Find largest eigenvaue of A using the power method.
% INPUT:
%    Afn - function handle for A
%    s - size of array that A acts on
%    Nit - number of iterations
% OUTPUT:
%    p - estimate of maximum eigenvalue of A
%
% % Example:
% N = 5;
% A = randn(N,N);
% Afn = @(x) A*x;
% Nit = 20;
% maxeig(Afn,[N 1],Nit)
%
% % compare with: 
% max(abs(eig(A)))

x = rand(s);
p = sqrt(sum(abs(x(:)).^2));
x = x/p;

if nargin < 3
    Nit = 10
end

for k = 1:Nit
    x = Afn(x);
    p = sqrt(sum(abs(x(:)).^2))
    x = x/p;
end


% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu


