% this code generates Figure 4 in Lutscher et al, Theoretical Ecology 2020
% parameter sweep for the simple pairformation model
% this program uses CalcSpeed.m for the speed calculation at each step

clear
tic
bmin = 1;  % minimum value of b
bmax = 3;  % maximum value of b 
db = 0.1;  % stepsize in b

nmin = 1;  % minimum value in n
nmax = 5;  % maximum value in n
dn = 0.1;  % stepsize in n

mu = 0;    % mu = 0 corresponds to no mutation

ll=0; kk=0; % counting variables

% calculating the speed in the absence of mutation
for b = bmin:db:bmax
  ll=ll+1;
  for n = nmax:-dn:nmin
    kk=kk+1;
    Bval(kk,ll)=b;
    Nval(kk,ll)=n;
    Speed(kk,ll) = CalcSpeed([b,n,mu]);
  end
  kk=0;
end
figure(1)
surf(Bval,Nval,max(Speed,0))
colorbar


mu = 0.1;  % mu = 0 corresponds to small mutation

ll=0; kk=0; % counting variables
% calculating the speed in the presence of mutation
for b = bmin:db:bmax
  ll=ll+1;
  for n = nmax:-dn:nmin
    kk=kk+1;
    Bval(kk,ll)=b;
    Nval(kk,ll)=n;
    Speed1(kk,ll) = CalcSpeed([b,n,mu]);
  end
  kk=0;
end
figure(2)
surf(Bval,Nval,max(Speed,0))
colorbar

% calculating the difference between the two

SpeedDiff1 = max(Speed1,0)-max(Speed,0);
figure(3)
surf(Bval,Nval,SpeedDiff1)
colorbar
toc