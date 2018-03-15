function [x1, x2, c1, c2, cost] = dualBP(x, A1, A1T, p1, A2, A2T, p2, lam1, lam2, mu1, mu2, Nit, plot_flag)

% [x1, x2, c1, c2] = dualBP(x, A1, A1T, p1, A2, A2T, p1, lam1, lam2, mu1, mu2, Nit);
%
% DUAL BASIS PURSUIT
% Signal decomposition using two tranforms (Parseval frames)
% Minimize:
% lam1 ||c1||_1 + lam2 ||c2||_2  subject to x = A1 c1 + A2 c2
%
% INPUT
%   x                : input signal signal
%   A1T, A1, A2T, A2 : functions handles for Parseval frames and inverses
%   p1, p2           : Parseval constants (A1*A1T = p1*I, A2*A2T = p2*I)
%   lam1, lam2       : regularization parameters
%   mu1, mu2         : ADMM parameters
%   Nit              : number of iterations
%
% OUTPUT
%   x1, x2 : components
%   c1, c2 : coefficients of components
%
% Use [x1,x2,c1,c2,cost] = dualBP(...) to return cost function history.
%
% Use [...] = dualBP(...,'plots') to plot progress of algorithm.

% Ivan Selesnick
% Polytechnic Institute of New York University
% May 2010
% Revised August 2011

% Reference:
% I. W. Selesnick. Sparse signal representations using the tunable Q-factor wavelet transform.
% Proc. SPIE 8138 (Wavelets and Sparsity XIV), August 2011. doi:10.1117/12.894280.

% The algorithm is a variant of SALSA (Afonso, Bioucas-Dias, Figueiredo,
% IEEE Trans Image Proc, 2010, p. 2345)


% By default do not compute cost function (to reduce computation)
if nargout > 4
    COST = true;
    cost = zeros(1,Nit);     % cost function
else
    COST = false;
    cost = [];
end

GOPLOTS = false;
if nargin == 13
    if strcmp(plot_flag,'plots')
        GOPLOTS = true;
    end
end

% Initialize:

c1 = A1T(x);
c2 = A2T(x);

d1 = A1T(zeros(size(x)));
d2 = A2T(zeros(size(x)));

T1 = lam1/(2*mu1);
T2 = lam2/(2*mu2);

N = length(x);
A = 1.1*max(abs(x));

C1 = (1/mu1)/(p1/mu1 + p2/mu2);
C2 = (1/mu2)/(p1/mu1 + p2/mu2);

for k = 1:Nit
    % fprintf('Iteration %d\n', k)
    
    u1 = soft(c1 + d1, T1) - d1;
    u2 = soft(c2 + d2, T2) - d2;
    
    c = x - A1(u1) - A2(u2);
    % c = 0.5 * c;
    
    d1 = C1*A1T(c);
    d2 = C2*A2T(c);
    
    c1 = d1 + u1;
    c2 = d2 + u2;
    
    if COST
        cost(k) = lam1*sum(abs(c1(:))) + lam2*sum(abs(c2(:)));
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

