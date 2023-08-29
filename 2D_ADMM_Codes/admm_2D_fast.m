function  im=admm_2D_fast(Fr,Fa,Y,e,alpha)


temp=1/sqrt(size(Fr,2))*fft(Y,size(Fr,2),1);
im=sqrt(size(Fa,2))*ifft(temp,size(Fa,2),2);


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


    temp1=sqrt(size(Fr,2))*ifft(B-U,size(Fr,2));
    temp1=temp1(1:size(Fr,1),:);
    temp2=1/sqrt(Na)*fft(temp1,[],2);
    temp2=temp2(:,1:size(Fa,1));
    temp3=temp2-Y;
    temp3=1/sqrt(Nr)*fft(temp3,Nr,1);
    temp3=sqrt(size(Fa,2))*ifft(temp3,size(Fa,2),2);

    im_k1= (B-U) -1/(delta+1)*temp3;

    B = shrink(im_k1 + U , alpha/delta);


    ERROR=norm(im_k1(:)-im(:))/norm(im(:));

    if ERROR<=e
        break
    end


    it=it+1;
    %      alphaIncrement = 1.02;

    U = U+ im_k1 - B;

    im=im_k1;

end

    function A = shrink(B,gamma)
        A = sign(B).*max(abs(B)-gamma,0);
    end
end