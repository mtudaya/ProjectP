function displaySTFT(c, fs, T, dblim)

% displaySTFT(c, fs, T, dblim)
%
% DISPLAY STFT (short-time Fourier transform)
%
% INPUT
%   c     : STFT coefficient (2D array)
%   fs    : sampling rate (samples/second)
%   T     : signal duration (seconds)
%   dblim : limits for colorbar in Db

[Nf, Nt] = size(c);

imagesc([0 T], [0 fs/2]/1e3, 20*log10(abs(c(1:Nf/2,:))), dblim)
axis xy
ylabel('Frequency (kHz)')
xlabel('Time (seconds)')

cm = colormap( 'gray' );
cm = cm(end:-1:1,:);
colormap(cm);

CB = colorbar;
set(CB, 'ytick', (-50:10:0))

shg

