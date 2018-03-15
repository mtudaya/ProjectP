function b = convt(h,g);
% b = convt(h,g);
% Transpose convolution: b = H' g.
% If h is complex, then H' is complex conjugate transpose

% Ivan Selesnick
% NYU-Poly
% selesi@poly.edu

Nh = length(h);
Ng = length(g);

b = conv(conj(h(Nh:-1:1)), g);
b = b(Nh:Ng);
