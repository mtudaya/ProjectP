% Signal separation with dual basis pursuit
%
% Ivan Selesnick

%% Start

clear
close all

addpath data

MyGraphPrefsON
% printme = @(filename) print('-deps', sprintf('figures/Example_dualBP_2_%s', filename) );
printme = @(filename) print('-dpdf', sprintf('figures/Example_dualBP_2_%s', filename) );
% printme = @(filename) disp(filename)

I = sqrt(-1);


%% Load data

[s,fs] = wavread('data/sp1.wav');
M = 20000;
y = s(1:M)';  

[y, fs] = wavread('data/arctic_a0001.wav' );
M = length(y);
y = y(:)';

m = (0:M-1)';

ymax = max(abs(y));

figure(1)
clf
subplot(2,1,1)
plot(m, y);
xlim([0 M])
ylim([-ymax ymax])
box off
mytitle('Signal')
xlabel(' ')

% printme('signal')

%% Define the two Parseval frames

if 0
% A1 : STFT
R1 = 128;
R1 = 256;
A1T = @(x) pSTFT(x, R1, R1);
A1 = @(X) ipSTFT(X, R1, M);
p1 = 1;

% A2 : STFT
R2 = 1024;
A2T = @(x) pSTFT(x, R2, R2);
A2 = @(X) ipSTFT(X, R2, M);
p2 = 1;
end

% A1 : STFT
R1 = 128;
% R1 = 256;
A1T = @(x) pSTFT2(x, R1, 4, 3, R1);
A1 = @(X) ipSTFT2(X, R1, 4, 3, M);
p1 = 1;

% A2 : STFT
R2 = 1024;
R2 = 512
A2T = @(x) pSTFT2(x, R2, 4, 3, R2);
A2 = @(X) ipSTFT2(X, R2, 4, 3, M);
p2 = 1;



% Check Parseval energy identity
E = sum(abs(y(:)).^2);
c1 = A1T(y);
c2 = A2T(y);
E1 = sum(abs(c1(:)).^2);
E2 = sum(abs(c2(:)).^2);
fprintf('Signal energy = %.3e\n', E)
fprintf('Transform 1 energy = %.3e\n', E1/p1)
fprintf('Transform 2 energy = %.3e\n', E2/p2)

% Check perfect reconstruction property
err = p1*y - A1(c1);
max(abs(err(:)))

err = p2*y - A2(c2);
max(abs(err(:)))



%% Show both STFT


dblim = [-50 0];

figure(1)
clf

subplot(2,1,1)
displaySTFT( c1, fs, M/fs, dblim)
title(sprintf('Original speech, STFT with %d-point window', R1))
% printme('STFT_y_R1')

subplot(2,1,2)
displaySTFT( c2, fs, M/fs, dblim)
title(sprintf('Original speech, STFT with %d-point window', R2))
% printme('STFT_y_R2')

orient tall
printme('STFT_y')

%% Peform signal separation
% Use the command 'dualBP' (dual basis pursuit)

theta = 0.5;
% theta = 0.45;
Nit = 200;
mu1 = 10.0;
mu2 = 10.0;


[y1,y2,c1,c2,costfn] = dualBP(y, A1, A1T, p1, A2, A2T, p2, theta, 1-theta, mu1, mu2, Nit);


%% Display cost function

figure(2)
clf
it1 = 1;
% plot(it1:Nit, costfn(it1:Nit), '.-')
plot(1:Nit, costfn)
xlim([0 Nit])
box off
title('Cost function')
xlabel('Iteration')

costfn(end)

% axis([0 Nit 640 740])

%% Display decomposition

y1 = A1(c1);
y2 = A2(c2);
res = y - y1 - y2;


figure(5)
clf
plot(m/fs, y, m/fs, y1-0.4, m/fs, y2-0.8)
legend('Speech waveform','Component 1','Component 2')
% xlim([0 M/fs])
% xlim([0.5 1])
xlim([2.0 2.5])
box off
xlabel('Time (seconds)')

orient landscape
printme('components')
orient portrait

%% Display coefficients

dblim = [-50 0];

g1 = A1T(y1);
g2 = A2T(y2);

g1 = g1/max(abs(g1(:)));
g2 = g2/max(abs(g2(:)));

figure(1)
clf
subplot(2,1,1)

% displaySTFT( g1, fs, M/fs, dblim)
displaySTFT( g1, fs, M/fs, [-40 0])
% title('Component 1')
title(sprintf('Component 1, STFT with %d-point window', R1))
% printme('STFT_y1_R1')

subplot(2,1,2)
displaySTFT( g2, fs, M/fs, dblim)
% title('Component 2')
title(sprintf('Component 2, STFT with %d-point window', R2))
% printme('STFT_y2_R2')

orient tall
printme('STFT_y1_y2')



%%

R = 512;
s1 = pSTFT(y1, R, R);
s2 = pSTFT(y2, R, R);

figure(1)
clf
subplot(2,1,1)
displaySTFT( s1, fs, M/fs, dblim)
title(sprintf('Component 1, STFT with %d-point window', R))
% printme('STFT_y1_R')

% figure(2)
% clf
subplot(2,1,2)
displaySTFT( s2, fs, M/fs, dblim)
title(sprintf('Component 2, STFT with %d-point window', R))
% printme('STFT_y2_R')

orient tall
printme('STFT_y1_y2_R')


%%

sound(y1, fs)

pause(.2)

sound(y2, fs)
