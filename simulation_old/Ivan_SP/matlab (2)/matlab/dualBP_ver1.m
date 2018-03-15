function [x1,x2,c1,c2,costfn] = dualBP_ver1(x, A1, A1T, A2, A2T, lam1, lam2, mu, Nit, plot_flag)
% Dual Basis Pursuit
% [x1,x2,c1,c2,costfn] = dualBP_ver1(x, A1T, A1, A2T, A2, lam1, lam2, mu, Nit)
% Signal decomposition using two tranforms (Parseval frames)
% This function minimizes the cost function:
%   lam1 * ||c1||_1 + lam2 * ||c2||_2  subject to x = A1(c1) + A2(c2)
%
% INPUT
%   x - input signal signal
%   A1T, A1, A2T, A2 - handles of two parseval frames and inverses
%   lam1, lam2 - parameters
%   mu - SALSA parameter
%   Nit - Number of iterations
% OUTPUT
%   x1, x2 - components
%   c1, c2 - coefficients of components
%
% Use [x1,x2,c1,c2,costfn] = dualBP(...) to return cost function.
%
% Use [...] = dualBP(...,'plots') to plot progress of algorithm.
%
%
% Ivan Selesnick
% Polytechnic Institute of New York University
% May 2010
% Revised August 2011
%
%
% References
%
% I. W. Selesnick. Sparse signal representations using the tunable Q-factor wavelet transform.
% Proc. SPIE 8138 (Wavelets and Sparsity XIV), August 2011. doi:10.1117/12.894280 (15 pages).


% By default do not compute cost function (to reduce computation)
if nargout > 4
    COST = true;
    costfn = zeros(1,Nit);     % cost function
else
    COST = false;
    costfn = [];
end

GOPLOTS = false;
if nargin == 10
    if strcmp(plot_flag,'plots')
        GOPLOTS = true;
    end
end

L = length(x);

% Initialize:

c1 = A1T(x);
c2 = A2T(x);

d1 = A1T(zeros(size(x)));
d2 = A2T(zeros(size(x)));

T1 = lam1/(2*mu);
T2 = lam2/(2*mu);

N = length(x);
A = 1.1*max(abs(x));

for k = 1:Nit
%     fprintf('Iteration %d\n', k)
    
    u1 = soft(c1 + d1, T1) - d1;
    u2 = soft(c2 + d2, T2) - d2;
    
    c = x - A1(u1) - A2(u2);
    c = 0.5 * c;
    
    d1 = A1T(c);
    d2 = A2T(c);
    
    c1 = d1 + u1;
    c2 = d2 + u2;
    
    if COST
        costfn(k) = lam1*sum(abs(c1(:))) + lam2*sum(abs(c2(:)));
    end
    
    if GOPLOTS
        x1 = A1(c1);
        x2 = A2(c2);
        res = x - x1 - x2;
        
        figure(gcf)
        clf
        subplot(3,1,1)
        plot(real(x1))
        xlim([0 N])
        ylim([-A A])
        title({sprintf('ITERATION %d',k),'COMPONENT 1'})
        box off
        subplot(3,1,2)
        plot(real(x2))
        xlim([0 N])
        ylim([-A A])
        box off
        title('COMPONENT 2')
        subplot(3,1,3)
        plot(real(res))
        xlim([0 N])
        ylim([-A A])
        title('RESIDUAL')
        box off
        drawnow
    end
    
end

x1 = A1(c1);
x2 = A2(c2);


