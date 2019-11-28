function y = partialFourier(idx,n,x,mode)
%PARTIALFOURIER  Partial Fourier operator
%
% Y = PARTIALFOURIER(IDX,N,X,MODE)

if mode==1
   z = fft(x) / sqrt(n);
   y = z(idx);
else
   z = zeros(n,1);
   z(idx) = x;
   y = ifft(z) * sqrt(n);
end
