clear
N = 200; Q = 100000; K = 100; J = 1; kbT = 1.5;
Tmin = 1.5; Tmax = 3; dT = 0.001;
a = ceil(rand(N,N)*2)*2 - 3;
for q = 1 : Q
    kbT = kbT + dT;
    if (kbT < Tmin) || (kbT > Tmax)
        dT = -dT;
    end
    kbT
    p = [1 1 1 exp(-4*J/kbT) exp(-8*J/kbT)];
    r0 = ceil(rand(N^2,2)*N); rn = mod(r0 - 2,N) + 1; rp = mod(r0,N) + 1;
    r = rand(N^2,1);
    for n = 1 : N^2
        if (r(n) < p(((a(rn(n,1),r0(n,2)) + a(rp(n,1),r0(n,2)) + a(r0(n,1),rn(n,2)) + a(r0(n,1),rp(n,2)))*a(r0(n,1),r0(n,2))/2 + 3)))
            a(r0(n,1),r0(n,2)) = -a(r0(n,1),r0(n,2));
        end
    end
    imagesc(a);
    axis equal off;
    drawnow;
end