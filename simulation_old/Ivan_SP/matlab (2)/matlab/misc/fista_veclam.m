function [x,J] = fista_veclam(y,H,Ht,lambda,alpha,Nit,x_init)
% [x, J] = fista_veclam(y,H,Ht,lambda,alpha,Nit)
%
% l1 regularization with 'fast iterated soft-thresholding algorithm'
% by Beck and Teboulle,
% Minimizes J(x) = ||y-H*x||_2^2 + ||lambda.*x||_1
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
    J(k) = sum(abs(Hx(:)-y(:)).^2) + sum(lambda(:).*abs(x(:)));
    %     J(k)

    Hg = H(g);
    x_old = x;
    x = soft(g + (Ht((y - Hg)))/alpha, T);

    t(k+1) = (1+sqrt(1+4*t(k)^2))/2;
    g = x + (t(k)-1)/t(k+1) * (x - x_old);

    if mod(k,10) == 0
        figure(101)
        clf
        CLim = [-70 0];
        C = max(abs(g(:)));
        ImageAD(fftshift(g)/C, CLim);
        title(sprintf('ITERATION %d. OBJECT FUN = %g',k,J(k)), 'FontSize',16)
        drawnow
    end

    PrintIterates = 0;
    if PrintIterates
        figure(101)
        clf
        CLim = [-70 0];
        C = max(abs(g(:)));
        ImageAD(fftshift(g)/C, CLim);
        title(sprintf('ITERATION %d. OBJECTIVE FUNCTION = %g',k,J(k)), 'FontSize',16)

        print('-djpeg',sprintf('FISTA_iterations/Iteration_%.3d',k))
    end

end
