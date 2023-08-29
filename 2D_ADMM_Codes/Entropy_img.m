function [ IE ] = Entropy_img( ISAR )
ISAR1=abs(ISAR).^2;
SumU = sum(sum(ISAR1));
II = (ISAR1/SumU);
Emat = II.*log(II);
IE = -(sum(sum(Emat)));
end

