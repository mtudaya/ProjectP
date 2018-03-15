% Relative amplitude distribution
% relative amplitude represents surface rougnness
% exponential distribution
%pdf(A)=1/mu exp (-A/mu);
function[A]=Relative_amp(A_mu,A_max,A_min)
A_mu=0.246;
pdf_A = makedist('Exponential','mu', A_mu);
A=0;

A=random(pdf_A);
if A>A_max || A< A_min
    
    while A > A_max || A< A_min
        A=random(pdf_A);
    end
end



return
