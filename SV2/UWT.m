% Compare DCTG2 and UWT, the Second Generation Discrete Curvelet 
% Transform and the Undecimated Wavelet Transform

% Read image;
img = fliplr(rot90(ReadImage('Barbara'),-1));
[N,J] = quadlength(img);
subplot(221);imagesc(img);axis image;axis off;colormap('gray')
     set(gca,'FontSize',14);
     title('(a)');

% Coarsest decomposition scale.
     coarsest = 2;

% Generate Symmlet 4 CMF Filter.
qmf = MakeONFilter('Symmlet',4);

% Compute UWT and DWT of noisy image.
uwc = FWT2_TI(img,coarsest,qmf);

% Plot sorted UWT coeffs.
swc   = sort(abs(uwc(:)),1,'descend');
subplot(222);plot(log2(1:(3*(J-coarsest)+1)*N*N),log2(swc),'--');
set(gca,'FontSize',14,'XTick',[0:2:log2((3*(J-coarsest)+1)*N*N)]);
xlabel('log_2(Index of retained DWT coeff)');ylabel('log_2(|coeff])');
title('(b)');

% Take DCTG2 transform.
C = fdct_wrapping(img,1,J-coarsest);

% DCTG2 implements a Parseval tight frame
% and the curvelets have l_2 norm  = 1/sqrt(redundancy).
Cvec = C{1}{1}(:);
nb=prod(size(C{1}{1}));
for j=2:length(C)
     for w=1:length(C{j})
     nb=nb+prod(size(C{j}{w}));
     Cvec = [Cvec;C{j}{w}(:)];
  end
end
     E = N/sqrt(nb);

% Plot sorted UWT coeffs.
  sC   = sort(abs(Cvec(:)/E),1,'descend');
subplot(222);hold on;plot(log2(1:nb),log2(sC));hold off;
axis([6 20 log2(sC(2^20)) log2(swc(2^6))]);
legend('UWT','DCTG2','Location','Best');
title('(b)');

% Keep only 2% of DCTG2 coefficients.
thd = sC(floor(nb*0.02));
Ct = C;
for j=2:length(C)
     for w = 1:length(C{j})
     Ct{j}{w} = C{j}{w}.* (abs(C{j}{w}) >= thd*E);
  end
end

% Take inverse curvelet transform
     imrdctg2 = real(ifdct_wrapping(Ct,1,J-coarsest));

subplot(223)
     imagesc(imrdctg2);axis image;axis off
set(gca,'FontSize',14);
title('(c)');
