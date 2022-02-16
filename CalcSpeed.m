function y = CalcSpeed(Params)

% Calculate the speed in the pairformation model
% with different variance of the kernel and different birth rate and
% carrying capacity

b = Params(1);
n = Params(2);
mu = Params(3);


% parameters
d1 = 1; % short dispersal parameter
d2 = 5; % long dispersal parameter
ctau = 1; % mate search rate
Tmax = 40; % number of generations

% set-up of the spatial grid
l = 200; np = 2^14;  dx = l/np;  x = linspace(-l/2,l/2-dx,np);    
xx = linspace(-l, l-dx, 2*np);   % double domain for padding
PAD = (abs(xx)<=l/2); % padding outside [-l/2, l/2]

% definition of the kernels and initial conditions
GAUSS1 = 1/sqrt(2*pi*d1)*exp(-xx.^2/(2*d1)); 
GAUSS2 = 1/sqrt(2*pi*d2)*exp(-xx.^2/(2*d2)); 
FG1 = fft(GAUSS1);
FG2 = fft(GAUSS2);
N10 = n/2*(xx<0); % total number of short distance dispersers
N20 = N10;

N1 = N10;
N2 = N20;
f = N1./(N1+N2+0.0000001);

for kk=1:Tmax
  N = N1+N2;
  P = (N.^2)./(N+1/ctau);
  f = N1./(N+0.0000001); 
  fnew = f+ mu*(1-2*f);
  N1new =  min(b*P,n).*fnew; 
  N2new = min(b*P,n).*(1-fnew);
 
  FF1 = fft(N1new); FF2 = fft(N2new);
  
  F1 = dx*real( fftshift( ifft( FF1.*FG1 ) ) ); 
  F2 = dx*real( fftshift( ifft( FF2.*FG2 ) ) ); 
   
  N1 = F1.*PAD.*(F1>0.000000001); 
  N2 = F2.*PAD.*(F2>0.000000001);
    
  if max(N1+N2)>0.05
  Front(kk) = dx*(max(find((N1+N2>0.01))))-l;
  else
  Front(kk)=-l;
  end
  
end

y = (Front(end)-Front(end-20))/20;

end
