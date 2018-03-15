function X = pSTFT(x,R,Nfft)

% X = pSTFT(x,R,Nfft)
%
% Parseval short-time Fourier Transform with 50% overlap.
% Each block is multiplied by a cosine window.
%
% INPUT
%   x    : signal (1D)
%   R    : block length (must be even)
%   Nfft : FFT  length (Nfft >= R)
%
% OUTPUT
%   X    : STFT coefficients (2D array)
%
% % EXAMPLE
%   [s,fs] = wavread('data/sp1.wav');
%   N = 20000;
%   x = s(1:N)';  
%   R = 512;
%   Nfft = 1024;
%   X = pSTFT(x,R,Nfft);
%   y = ipSTFT(X,R,N);
%   max(abs(x - y))       % verify perfect reconstruction
%   sum(abs(x(:)).^2)     % signal energy
%   sum(abs(X(:)).^2)     % STFT energy (equal to signal energy)
%   displaySTFT(X,fs,N/fs,[-50 0]);


x = x(:).';                             % Ensure x is row vector
n = (1:R) - 0.5;
win  = sin(pi*n/R);                     % cosine window
M = 2;
NC = sqrt(sum(win.^2) * M * Nfft/R);    % normalization constant
x = [zeros(1,R) x zeros(1,R)];          % to deal with first and last block
X = buffer(x,  R, R/2, 'nodelay');
X = bsxfun(@times, win', X);
X = fft(X, Nfft);                       % FFT applied to each column of X
X = X/NC;

