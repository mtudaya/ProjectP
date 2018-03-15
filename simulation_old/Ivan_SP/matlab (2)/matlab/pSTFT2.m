function X = pSTFT2(x,R,M,K,Nfft)

% X = pSTFT2(x,R,M,K,Nfft)
%
% Parseval short-time Fourier Transform with hop fraction 1/M.
% Intended use: STFT with more than 50% overlap.
% [overlap fraction = (1-1/M)].
% - M-times overcomplete
% - window = sin(pi*n/R).^K
%
% INPUT
%   x    : signal (1D)
%   R    : block length
%   M    : over-sampling rate
%   K    : window shape parameter
%   Nfft : FFT length (Nfft >= R)
% NOTES
%   R should be an integer multiple of M.
%   K should be less than M (1 <= K < M).
%
% OUTPUT
%   X    : STFT coefficients (2D array)
%
% % EXAMPLE
%   [s,fs] = wavread('data/sp1.wav');
%   N = 20000;
%   x = s(1:N)';
%   R = 501;
%   M = 3;
%   K = 2;
%   Nfft = 512;
%   X = pSTFT2(x,R,M,K,Nfft);
%   y = ipSTFT2(X,R,M,K,N);
%   max(abs(x - y))       % verify perfect reconstruction
%   sum(abs(x(:)).^2)     % signal energy
%   sum(abs(X(:)).^2)     % STFT energy (equal to signal energy)
%   displaySTFT(X,fs,N/fs,[-50 0]);

x = x(:).';                             % ensure x is row vector
n = (1:R) - 0.5;
win  = sin(pi*n/R).^K;                  % cosine window
NC = sqrt(sum(win.^2) * M * Nfft/R);    % normalization constant
x = [zeros(1,R) x zeros(1,R)];          % to deal with first and last block
X = buffer(x,  R, R*(M-1)/M, 'nodelay');
X = bsxfun(@times, win', X);
X = fft(X, Nfft)/NC;                    % FFT applied to each column of X

