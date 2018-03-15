% Denoising using DCTG2, Second Generation Discrete Curvelet Transform 

% Read image;
img = fliplr(rot90(ReadImage('Barbara'),-1));
[N,J] = quadlength(img);
% Add AWGN.
imn = img + 20*randn(size(img));

% Coarsest decomposition scale.
coarsest = 2;

% Generate Symmlet 4 CMF Filter.
qmf = MakeONFilter('Symmlet',4);

% Compute UWT of noisy image.
uwc = FWT2_TI(imn,coarsest,qmf);

% Estimate noise std from diagonal band at finest scale.
hh = uwc(2*N+1:3*N,:);
sigma = MAD(hh(:));

% Hard-threshold UWT coeffs: finest scale at 4*sigma and other scales at 3*sigma.
uwc(1:3*N,:)       = HardThresh(uwc(1:3*N,:),4*sigma);       % Finest.
uwc(3*N+1:end-N,:) = HardThresh(uwc(3*N+1:end-N,:),3*sigma); % Other scales.

% Reconstruct.
imdenuwt = IWT2_TI(uwc,coarsest,qmf);

% Take DCTG2 transform.
C = fdct_wrapping(imn,1,J-coarsest);

% DCTG2 implements a Parseval tight frame
% and the curvelets have l_2 norm  = 1/sqrt(redundancy).
nb=prod(size(C{1}{1}));
for j=2:length(C)
     for w=1:length(C{j})
     nb=nb+prod(size(C{j}{w}));
  end
end
E = N/sqrt(nb);

% Hard-threshold UWT coeffs: finest scale at 4*sigma and other scales at 3*sigma.
  Ct = C;
for j=2:length(C)
     thresh = 3*sigma + sigma*(j == length(C));
     for w=1:length(C{j})
     Ct{j}{w} = C{j}{w}.* (abs(C{j}{w}) > thresh*E);
  end
end

% Take inverse curvelet transform
     imdendctg2 = real(ifdct_wrapping(Ct,1,J-coarsest));

% Display.
  subplot(221)
     imagesc(img);axis image;axis off;colormap('gray')
     set(gca,'FontSize',14);
     title('(a)');
     subplot(222)
     imagesc(imn);axis image;axis off
     set(gca,'FontSize',14);
title('(b)');
subplot(223)
     imagesc(imdendctg2);axis image;axis off
set(gca,'FontSize',14);
title('(c)');
subplot(224)
     imagesc(imdenuwt);axis image;axis off
set(gca,'FontSize',14);
title('(d)');
