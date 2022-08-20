clear
Tc = 2/log(1+sqrt(2));
N = 200; K = 7000; J = 1;
T = (0.5:0.1:1.5)'*Tc;
%T = flip(T);  % to move from higher temperature to lower temperature

M(size(T)) = 0;
a = ceil(rand(N,N)*1.5)*2 - 3;
for t = 1 : size(T,1)
    p = [1 1 1 exp(-4*J/T(t)) exp(-8*J/T(t))];
    s(1:K,1) = 0;
    for k = 1 : K
        r0 = ceil(rand(N^2,2)*N); rn = mod(r0 - 2,N) + 1; rp = mod(r0,N) + 1;
        r = rand(N^2,1);
        for n = 1 : N^2
            if (r(n) < p(((a(rn(n,1),r0(n,2)) + a(rp(n,1),r0(n,2)) + a(r0(n,1),rn(n,2)) + a(r0(n,1),rp(n,2)))*a(r0(n,1),r0(n,2))/2 + 3)))
                a(r0(n,1),r0(n,2)) = -a(r0(n,1),r0(n,2));
            end
        end
        s(k) = sum(sum(a))/N^2;
    end
    %M(t) = sum(s(K-99:K))/100; % calculate magnetization using last hundred steps
    M(t) = sum(s(1:100))/100; % calculate magnetization using last hundred steps
    imagesc(a);
    axis equal off;
    drawnow;
    [t M(t)]
end
figure;
plot(T/Tc,abs(M),'o');
hold on;
Tf = (0.5:0.001:1)*Tc';
plot(Tf/Tc,real((1 - sinh(2*J*Tf.^-1).^-4).^(1/8)));
set(gca,'FontSize',20);
xlabel('T/Tc');
ylabel('Magnetization');