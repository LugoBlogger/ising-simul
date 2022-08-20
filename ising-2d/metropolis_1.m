clear
N = 100; Q = 10000; J = 1; kbT = 2;
a = ceil(rand(N,N)*2)*2 - 3;
p = [1 1 1 exp(-4*J/kbT) exp(-8*J/kbT)];
for q = 1 : Q
    r0 = ceil(rand(N^2,2)*N); rn = mod(r0 - 2,N) + 1; rp = mod(r0,N) + 1;
    r = rand(N^2,1);
    for n = 1 : N^2
        if (r(n) < p( ...
            (a(rn(n,1),r0(n,2)) ...
                + a(rp(n,1),r0(n,2)) ...
                + a(r0(n,1),rn(n,2)) ...
                + a(r0(n,1),rp(n,2))) * a(r0(n,1),r0(n,2))/2 + 3))
            a(r0(n,1),r0(n,2)) = -a(r0(n,1),r0(n,2)); % flip the spin
        end
    end
    imagesc(a);
    axis equal off;
    drawnow;
end