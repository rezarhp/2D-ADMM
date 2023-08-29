function s=GP_SOONE(A,B, x, sigma_min, sigma_decrease_factor,L0,L1, L, A_pinv,B_pinv, true_s)

if nargin < 6
    sigma_decrease_factor = 0.5;
    B_pinv=pinv(B);

    A_pinv = pinv(A);
    L0 = 25;
    L1=10;
    L = 3;
    ShowProgress = false;
elseif nargin == 6
    B_pinv=pinv(B);

    A_pinv = pinv(A);
    L0 = 25;
    L1=10;
    L = 3;
    ShowProgress = false;
elseif nargin == 8
    B_pinv=pinv(B);
    A_pinv = pinv(A);
    ShowProgress = false;
elseif nargin == 9
    B_pinv=pinv(B);
    ShowProgress = false;
elseif nargin == 10
    ShowProgress = false;
elseif nargin == 11
    ShowProgress = true;
else
    error('Error in calling GP_SOONE function');
end

% Initialization
%s = A\x;
s = A_pinv*x*B_pinv.';
sigma = 100*max(abs(s(:)));
% Main Loop
J=round(-log(sigma/sigma_min)/log(sigma_decrease_factor))+1;
for j=1:J
    beta=(J-j/2+1)/J;
    for i=1:L
        gamma=(L-i/2+1)/L;
        mu_0=beta*gamma*min(max(abs(s(:)))/L0,sigma/L1);
        delta = OurDelta(s,sigma);
        s = s - mu_0*delta;
        s = s - A_pinv*(A*s*B.'-x)*B_pinv.';   % Projection

    end

    if ShowProgress
        fprintf('     sigma=%f, SNR=%f\n',sigma,estimate_SNR(s,true_s))
    end

    sigma = sigma * sigma_decrease_factor;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta=OurDelta(s,sigma)

delta = s./(abs(s)).*exp(-abs(s)/sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SNR=estimate_SNR(estim_s,true_s)

err = true_s - estim_s;
SNR = 10*log10(sum(abs(true_s).^2)/sum(abs(err).^2));