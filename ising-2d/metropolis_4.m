% metropolis_4.m
%
%   This script will compute 3D-ising 

clear
N = 100; 
Q = 10000; 
J = 1; 
kbT = 2;

a = (rand(N,N,N) < 0.5)*2 - 1;
p = [1 1 1 1 exp(-4*J/kbT) exp(-8*J/kbT) exp(-12*J/kbT)];

for q = 1:Q
  r0 = randi([1 N],N^3,3); rn = mod(r0 - 2,N) + 1; rp = mod(r0,N) + 1;
  r = rand(N^3,1);
  for n = 1 : N^3
      if (r(n) < p(((a(rn(n,1),r0(n,2),r0(n,3)) + a(rp(n,1),r0(n,2),r0(n,3)) + ...
                     a(r0(n,1),rn(n,2),r0(n,3)) + a(r0(n,1),rp(n,2),r0(n,3)) + ...
                     a(r0(n,1),r0(n,2),rn(n,3)) + a(r0(n,1),r0(n,2),rp(n,3)))*...
                     a(r0(n,1),r0(n,2),r0(n,3))/2 + 4)))
          a(r0(n,1),r0(n,2),r0(n,3)) = -a(r0(n,1),r0(n,2),r0(n,3)); % flip the spin
      end
  end
  imagesc(a(:,:,mod(q,N)+1));
  axis equal off;
  drawnow;
end