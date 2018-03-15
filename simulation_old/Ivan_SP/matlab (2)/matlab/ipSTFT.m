function y = ipSTFT(X,R,N)
% y = ipSTFT(X,R,N)
%
% Inverse Short-time Fourier Transform with 50% overlap.
% This is the inverse of 'pSTFT'
%
% INPUT
%   X : STFT produced by 'pSTFT'
%   R : block length
%   N : length of signal to be reconstructed
%
% OUTPUT
%   y : signal

[Nfft, Nc] = size(X);                   % get sizes

n = (1:R) - 0.5;
win  = sin(pi*n/R);                     % cosine window

M = 2;
NC = sqrt(sum(win.^2) * M * Nfft/R);    % normalization constant

Y = ifft(X);                            % inverse FFT of each column of X
Y = Y(1:R,:);                           % truncate down to block length
Y = bsxfun(@times, Y, win');

y = zeros(1,R/2*(Nc+1));
i = 0;
for k = 1:Nc
    y(i + (1:R)) = y(i + (1:R)) + Y(:,k).';
    i = i + R/2;
end

y = NC * y(R+(1:N));

