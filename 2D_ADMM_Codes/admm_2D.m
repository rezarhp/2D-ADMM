function  im=admm_2D(Fr,Fa,Y,e,alpha)

im=Fr'*Y*conj(Fa);

Na = size(im,2);
Nr = size(im,1);

B = zeros(Nr,Na);
U = zeros(Nr,Na);


maxIter = 40;

delta = 1;

it = 0;


cgtol=1e-4;
% alpha=0.1;


ERROR = 1;
while ERROR > cgtol  && it < maxIter


    im_k1= (B-U) -1/(delta+1)*Fr'*(Fr*(B-U)*Fa.'-Y)*conj(Fa);

    B = shrink(im_k1 + U , alpha/delta);


    ERROR=norm(im_k1(:)-im(:))/norm(im(:));

    if ERROR<=e
        break
    end


    it=it+1;

    U = U+ im_k1 - B;

    im=im_k1;

end

    function A = shrink(B,gamma)
        A = sign(B).*max(abs(B)-gamma,0);
    end
end