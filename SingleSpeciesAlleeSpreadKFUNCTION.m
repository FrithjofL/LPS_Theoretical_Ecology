function y = SingleSpeciesAlleeSpreadKFUNCTION(P)

% Lutscher, Popovic, Shaw, Theoretical Ecology
% this function is called by AllSpreadKTypesMortalityRUN.m
% this function runs a single type in isolation with a half line initial condition and returns the speed.

% parameters
b = P(1);  % offspring per female
n = P(2); % maximum local density
tau = P(3); % mate search rate
dmin = P(4); % minimal dispersal parameter 
dmax = P(5); % maximal dispersal parameter
m = P(8); % mortality factor
K = P(9); % number of types
Tmax = 100; % number of generations

% set-up of the spatial grid
l = 200; np = 2^15;  dx = l/np;  x = linspace(-l/2,l/2-dx,np);    
xx = linspace(-l, l-dx, 2*np);   % double domain for padding
PAD = (abs(xx)<=l/2); % padding outside [-l/2, l/2]
Dispersal = linspace(dmin,dmax,K);

dcount = 0;
for d = Dispersal
 dcount = dcount+1;
  % definition of the kernel and initial condition
  GAUSS1 = exp(-d*m)/sqrt(2*pi*d)*exp(-xx.^2/(2*d)); 
  FG1 = fft(GAUSS1);
  N = n*(xx<0);
  
  %figure(10)
  %plot(xx,N); axis([-l/2,l/2, 0, n+1]);
  
  for kk=1:Tmax
    F = min( b*((N).^2)./((N)+1/tau), n);
    FF = fft(F);   
    F = dx*real( fftshift( ifft( FF.*FG1 ) ) ); 
    N = F.*PAD.*(F>0.000000001);   
    if max(N)>0.1
      Front(kk) = dx*(max(find((N>0.01))))-l;
    else
      Front(kk)=0;
    end
    figure(10)
      if mod(kk,10)==0
      hold on 
      plot(xx,F,'r--')
      pause(0.01)
      hold off;
      end
   end
   % figure(11)
   % plot(1:Tmax,Front)
    speed = (Front(end)-Front(end-20))/20;
    y(dcount) = speed
  end
