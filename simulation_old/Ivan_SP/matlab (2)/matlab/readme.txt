
Matlab software to accompany the tutorial: 'Introduction to Sparsity in Signal Processing'

http://cnx.org/content/m43545/latest/
http://eeweb.poly.edu/iselesni/teaching/lecture_notes/sparsity_intro/index.html

Examples:
	Example_BP		basis pursuit example (sparse Fourier coefficients)
	Example_BPD		basis pursuit denoising example (speech denoising)
	Example_deconv		deconvolution using BPD
	Example_missing		missing data estimation using BP
	Example_dualBP_1	signal component separation (spikes + sinusoids)
	Example_dualBP_2	signal component separation (short + long STFT windows)

Matlab programs:
	A			oversampled DFT
	AT			conjugate transpose of 'A'
	bp_salsa		basis pursuit (BP) using algorithm SALSA
	bpd_salsa		basis pursuit denoising (BPD) using algorithm SALSA
	bpd_salsa_sparsemtx	implementation of BPD with sparse matrix A
	bp_missing		estimate missing data using BP
	dualBP			dual basis pursuit
	soft			soft thresholding
	pSTFT			Parseval STFT, 50% overlapping
	pSTFT2			Parseval STFT, flexible overlap factor
	ipSTFT			inverse of 'pSTFT'
	ipSTFT2			inverse of 'pSTFT2'
	displaySTFT		display STFT coefficients

Utility functions:
	MyGraphPrefsON		modify Matlab default graphing preferenes
	MyGraphPrefsOFF		set graphing preference back to Matlab default
	mytitle			variation on Matlab title function

Folders:
	data			data files for examples
	figures			figures produced by the examples for the tutorial
	html			html files produced by Matlab 'publish' function

Ivan Selesnick
Polytechnic Institute of New York University
selesi@poly.edu

April 2012

Support from NSF under Grant CCF-1018020 is gratefully acknowledged.


