function s=SL0_2D(A,B, X, sigma_min, sigma_decrease_factor, mu_0, L, A_pinv, true_s)


if nargin < 5
    sigma_decrease_factor = 0.5;
    A_pinv = pinv(A);
    mu_0 = 2;
    L = 3;
    ShowProgress = false;
elseif nargin == 5
    A_pinv = pinv(A);
    mu_0 = 2;
    L = 3;
    ShowProgress = false;
elseif nargin == 6
    A_pinv = pinv(A);
    L = 3;
    ShowProgress = false;
elseif nargin == 7
    A_pinv = pinv(A);
    ShowProgress = false;
elseif nargin == 8
    ShowProgress = false;
elseif nargin == 9
    ShowProgress = true;
else
    error('Error in calling SL0 function');
end

B_pinv=pinv(B);

% Initialization

s = A_pinv*X*B_pinv.';

sigma = 2*max(abs(s(:)));

% Main Loop
ii=0;
while sigma>sigma_min
    ii=ii+1;
    for i=1:L
        delta = OurDelta(s,sigma);
        s = s - mu_0*delta;
        s = s - A_pinv*(A*s*B.'-X)*B_pinv.';   % Projection
    end

    if ShowProgress
        fprintf('     sigma=%f, SNR=%f\n',sigma,estimate_SNR(s,true_s))
    end

    sigma = sigma * sigma_decrease_factor;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta=OurDelta(s,sigma)

delta = s.*exp(-abs(s).^2/sigma^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SNR=estimate_SNR(estim_s,true_s)

err = true_s - estim_s;
SNR = 10*log10(sum(abs(true_s).^2)/sum(abs(err).^2));