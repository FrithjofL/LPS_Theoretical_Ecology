function y = AllSpreadKTypesMortalityFUNCTION(P)

% Lutscher, Popovic, Shaw, Theoretical Ecology
% this function is called by AllSpreadKTypesMortalityRUN.m
% this function runs the model with K types


% parameters
b = P(1);   % offspring per female
n = P(2);  % maximum local density
tau = P(3);  % mate search rate
dmin = P(4);  % minimal dispersal parameter
dmax = P(5);  %maximal dispersal parameter
mu = P(6);   % mutation probability for simulation run
mui = P(7);  % mutation probability for  initial conditions
m = P(8); % mortality factor: survival is exp(-d(i)*m) 
K = P(9); % number of types
Tmax = 500; 
if mu==0
  Tmax = 1200; % number of generations
end

  
% set-up of the spatial grid
l = 500; np = 2^16;  dx = l/np;  x = linspace(-l/2,l/2-dx,np);    
xx = linspace(-l, l-dx, 2*np);   % double domain for padding
PAD = (abs(xx)<=l/2); % padding outside [-l/2, l/2]

% initial conditions
N0 = (xx<-150);  % half line 
%01 = (xx<-100);  % half line with advantage
% Dispersal vector 
D = linspace(dmin,dmax,K);
% Survival vector
DIAG = diag(exp(-D*m));
% mutation matrix_type for initial conditions only
MI = (1-mui)*eye(K);
for i=1:K-1
  MI(i,i+1)=mui/2;
  MI(i+1,i)=mui/2;
end
MI(2,1)=mui;
MI(K-1,K)=mui;

% mutation matrix_type for simulations only
M = (1-mu)*eye(K);
for i=1:K-1
  M(i,i+1)=mu/2;
  M(i+1,i)=mu/2;
end
M(2,1)=mu;
M(K-1,K)=mu;

% Initial conditions as a steady state
NN = n/K*ones(1,K)'; % initial assumption uniform distribution
for kkk=1:2000;
NNtot = sum(NN);
FF = NN/NNtot;
FFnew = MI*FF;
NNnew = min(b*NNtot^2/(NNtot+tau),n);
NN = NNnew*FFnew;
NN = DIAG*NN;
end


% initial conditions for joint spread simulations
for i=1:K
  N(i,:)=NN(i)*N0;
end

col = ['r','b','k','g','c','m','y','r','b','k','g','c','m','y'];
figure(1)
for i = 1:K
subplot(K,1,i), plot(xx,N(i,:),col(i)); axis([-l/2,l/2, 0, n+1]);
end

Nplot = zeros(size(xx));
figure(2)
for i=1:K
    Nplot = Nplot+N(i,:);
 %subplot(K+1,1, K+1), 
 plot(xx,Nplot,col(i)); axis([-l/2,l/2, 0, n+1]);
 hold on
end
hold off


% dispersal kernels and their FFTs
FG = zeros(K,length(xx));
for i=1:K
  GAUSS = exp(-D(i)*m)/sqrt(2*pi*D(i))*exp(-xx.^2/(2*D(i))); 
  FG(i,:) = fft(GAUSS);
  clear GAUSS
end

F = zeros(K,length(xx));
for tt=1:Tmax
  Ntot = sum(N);
  Nnew = min( b*((Ntot).^2)./((Ntot)+1/tau), n);
  for i=1:K
    F(i,:) = N(i,:)./(Ntot+0.0000001);
  end
    if max(Ntot)>0.05
      Front(tt) = dx*(max(find((Ntot>0.1))))-l;
    else
      Front(tt)=-l;
    end
  Fnew = M*F;
  N = Nnew.*Fnew;
  for i=1:K
    NN = dx*real( fftshift( ifft( fft(N(i,:)).*FG(i,:) ) ) ); 
    Nnext(i,:) = NN.*(NN>0.000000001).*PAD; 
  end
  N = Nnext;
    
  if mod(tt,50)==0
    figure(1)
    for i = 1:K
    hold on
    subplot(K,1,i), plot(xx,N(i,:),col(i)); axis([-l/2,l/2, 0, n+1]);
    hold off
    end
       
    figure(2)
    Nplot = zeros(size(xx));
    hold on
    for i=1:K
      Nplot = Nplot+N(i,:);
      %subplot(K+1,1,K+1), 
      plot(xx,Nplot,col(i)); axis([-l/2,l/2, 0, n+1]);
      hold on
    end
    
   
    % % for plot of two types stacked waves only
    %if mod(tt,100)==0
    %figure(10)
    %subplot(2,1,1)
    %hold on
    %plot(xx,N(1,:),'b--');
    %plot(xx,N(2,:),'b-');
    %axis([-l/2,l/2, 0, n+1]);
    %hold off
    %set(gca, "linewidth", 1, "fontsize", 12)
    %xlabel('space',"fontsize",16)
    %ylabel('density',"fontsize",16)
    %end
   
   
  hold off
  end


end

figure(1)
for i = 1:K
      subplot(K,1,i)
      hold off
end

    figure(3)
    hold on
    plot(xx,N(2,:)./(N(1,:)+N(2,:)+0.00001)); axis([-l/2,l/2, 0, 1])
    hold off 

%figure(11)
%plot(1:Tmax,Front)
y = (Front(end)-Front(end-20))/20;





