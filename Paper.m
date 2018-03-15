close all; clear; clc;
M = 6; % Number of antennas
B = 20e+6; % bandwidth of the receiver
Fs = 40e+6; % sampling rate (2*B)
Fc = 2.4e+9; % carrier frequency
c = 3e+8; % speed of EM
lambda = c/Fc; % wavelength
Ts = 1/Fs; % sampling interval
L = 6; % number of reflectrots
N =20; %length of the channel (-N to N)
%alpha = sort (0 + (1-0).*rand(L,1),'descend') % r = a + (b-a).*rand(N,1)
alpha=[0.8 0.7 0.5 0.3 0.1 0.05]; % reflection coefficients of the reflectors
d = lambda/2; % antenna spacing of the antenna array
theta = [80 -60 40 -25 10 -45]; % angle of arrival of the reflectors
v = [10 10 10 10 10 10];
T = [(0.25*Ts) (0.75*Ts) (Ts+0.5*Ts) (2*Ts+0.25*Ts) (3*Ts+0.25*Ts) (5*Ts+0.25*Ts)];
%theta = randi([-90 90],1,L)
%gamma = (2*pi/lambda)*d*(m-1)*(sin(pi*(theta(m))/180));
delay_min=0;
delay_max=10*Ts;
ND = 20; %resolution of the dictionary
p = (delay_max-delay_min)/ND; %size of the dictionary
 gamma = ones(M,L);
 for m = 1:M
        for l = 1:L
        gamma(m,l) = (2*pi/lambda)*d*(m-1)*(sin(pi*(theta(l))/180));
        end
 end


a = zeros(1,2*N+1);
b = zeros(M,2*N+1);
channel_RX_ANT = complex(a,0);
channel_Total = complex(b,0);
n = -N:1:N;
for m = 1:1:M
    
    a = zeros(1,2*N+1);
    channel_RX_ANT = complex(a,0);
    %disp(m);  
    for k = 1:1:L
        %disp(k);
        %disp(channel);
        channel_RX_ANT = channel_RX_ANT + (alpha(k)* exp(1i*(v(k)*pi/180)+gamma(m,k))*sinc(B*((n*Ts)-T(k)-(gamma(m,k)/(2*pi*Fc)))));
    end
    disp(channel_RX_ANT);
    channel_Total(m,:)=channel_RX_ANT;
    
end
%plot (n,channel_RX_ANT);hold on;
%plot (n,channel_Total);

% optimization section
% consider an one antenna

y = channel_Total(1,:);
b= transpose(y);
% construction of dictionaty D (nxP)
D = zeros(2*N+1,ND);
X1 = zeros(ND,1);
X = complex(X1,0);

% filling the D matrix

for index =-N:1:N
    
    for s = 1:1:ND
        disp(s)
        D((index+N+1),s)= sinc((B*((index*Ts)-(s*p))));
    end
end
 
%optimization part

% Generate an initial point.
x0       = 0.*ones(ND, 1);

%% Test the BP problem.
tolx = 1e-6;

% Set the parameters.
param.MaxIters      = 3000;
param.Verbosity     = 2;
param.RelTolX       = tolx;
param.saveHistMode  = 1;
param.Algorithm     = 1;  % You should choose 1 or 3 ...
param.InnerMaxIters = 2;  % This option can be set to 5 -> 10 ...
param.adaptStepSize = 0;  

% Call the solver - By default ...
if 1
    [x1, out1] = decoptSolver('BP', D, b, param, 'x0', x0);
else
    % User-define proximal-functions.
    proxOpers{1} = @(x, gamma, varargin)(proxL1norm(x, gamma));
    proxOpers{2} = @(x, gamma, varargin)(projL2norm(x, 1e-12));

    proxOpers{3} = @(x, varargin)(norm(x(:), 1));
    proxOpers{4} = @(x, varargin)(0);

    % Call the solver with user-define prox-functions.
    [x1, out1] = decoptSolver('UserDef', D, y, param, 'x0', x0, 'Prox', proxOpers);
end

% Evaluate the objective values and feasibility gap.
fx1   = norm(x1, 1);
feas1 = norm(D*x1 - y, 2); 

%% Printing ...
fprintf('******************** THE FINAL RESULTS ************************\n');
fprintf('+ DECOM: BP-problem: [f(x), |A*x-b|/|b|] = [%3.7f, %5.7f]\n', fx1, feas1/norm(b));
fprintf('+ DECOM: Iterations: %4d, Time(s) = %3.4f\n', out1.iter, out1.total_time);
fprintf('+ DECOM: Number of Ax and ATy are %4d and %4d\n', out1.cntA, out1.cntAt);
%fprintf('+ DECOM: Reconvery error: %4.7f\n', norm(x1 - x_org)/max(norm(x_org), 1));


%% Plot the figures.
%if isPlotFigure 
    
    % Plot the solution.
    %figure(1); title('The solutions');
    %if isreal(x_org), stairs(x_org, 'g:*');  else stairs(abs(x_org), 'g:*'); end
    %hold on;
    %if isreal(x1), stairs(x1, 'r--o');  else stairs(abs(x1), 'r--o'); end
    %hold off;
    %shg;
%end

%% END OF THE TEST.