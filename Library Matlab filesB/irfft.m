function x = irfft(X,nfft)
%IRFFT inverse (scaled) real fft (complex to real).
%  irfft(X) is the inverse discrete Fourier transform of vector X which
%     is known to be symmetric so that the result has to be real.
%  irfft(X,N) is the corresponding N-point inverse transform.
% X is assumed scaled as in rfft().
% See also FAST(), FSST() by G. Kubin 08-28-92, Jont Allen 9-25-92
% C.A.Shera 6-19-02

  [k,l]=size(X);
  if (k==1)
    nf=l;
    % zero out Imag part in DC and Fmax
    X(:,1)=real(X(:,1));
    X(:,nf)=real(X(:,nf));
    X=[X,conj(X(:,nf-1:-1:2))];
  else
    nf=k;
    % zero out Imag part in DC and Fmax
    X(1,:)=real(X(1,:));
    X(nf,:)=real(X(nf,:));
    X=[X;conj(X(nf-1:-1:2,:))];
  end

  if (nargin == 2)
    x = ifft(X,nfft);
  else
    x = ifft(X);
  end
  x = real(x) + imag(x);

  [m,n] = size(x);
  if m == 1
    m = n;
  end
  x = (m/2) * x;
