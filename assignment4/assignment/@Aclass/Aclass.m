function A=Aclass(PSF,center,bc,n,m)

if strcmp(bc,'periodic')
    radius=0;
else
    radius=ceil(max(size(PSF)/2));
end

PSFp=padarray(PSF,[n+2*radius,m+2*radius]-size(PSF),0,'post');
eigA=fft2(circshift(PSFp,1-center));

A.radius=radius;
A.eigA=eigA;
A.bc=bc;
A.n=n;
A.m=m;
A.PSF=padarray(PSF,[n,m]-size(PSF),0,'post');
A.center=center;
A=class(A,'Aclass');
end