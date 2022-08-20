% metropolis_1_v2.m
%
%  This script will compute magnetization for given temperature.
%  k_B is set to 1.
%  Boundary condition is periodic.

%  Logs
% 
%  [22/05/29]
%  - improve the computation of magnetization 
%
%  [22/06/19]
%  - add observable of specific heat
%
%  [22/07/09]
%  - revision to the computation of specific heat
%    Computation of specific heat can be done by computing the standard
%    deviation of energy.
%    In the previous version, I made a mistake by calculating specific heat
%    for a single Monte Carlo sweep. It should be calculated by taking
%    several sample after equilibration and then take the standard
%    deviation from those samples. 
%    Therefore the possible observable that we can compute is energy for
%    each Monte Carlo sweep

clear;

%% -- user-defined input
N = 50; 
L = 2500;   % number of samples after one MCS (Monte Carlo Sweep)
J = 1; 
kb = 1;
T_i = 3;   % Tc = 2.269

kbT = kb*T_i;   % we set k_B to 1

% initial configuration of spins are random
x_config_spin = ceil(rand(N,N)*2)*2 - 3;    % [0, 2] -> {1, 2} -> {2, 4} -> {-1, 1};

% probability that guides each flip is allowed or not
p = [1, ...              % Delta H = -8J < 0
     1, ...              % Delta H = -4J < 0
     1, ...              % Delta H = 0 -> p(0) = exp(0) = 1
     exp(-4*J/kbT), ...  % Delta H = 4J  [allowed with probability]
     exp(-8*J/kbT)];     % Delta H = 8J  [allowed with probability]

%% -- run Monte Carlo simulation for 2D Ising model
tic;

% -- ising simulation figure handle
fig_ising = figure('Units', 'centimeters', 'Position', [3, 3, 12, 12], 'Color', 'w');
  ax_ising = axes('Parent', fig_ising);
   
% -- magnetization curve figure handle
fig_magnetization = figure('Units', 'centimeters', 'Position', [16, 3, 20, 12], 'Color', 'w');
  ax_magnetization = axes('Parent', fig_magnetization); 
  exact_solutionHandler = plot(ax_magnetization, ...
    [0, L], ones(1, 2)*get_magnet_exact(kbT, J), ...
    'LineWidth', 1.5);
  hold(ax_magnetization, 'on');
  
  ylim(ax_magnetization, [-0.1, 1.1]);
  xlabel(ax_magnetization, 'MCS');
  ylabel(ax_magnetization, 'total |S| / N');

% -- energy curve figure handle
fig_energy = figure('Units', 'centimeters', 'Position', [17, 4, 20, 12], 'Color', 'w');
  ax_energy = axes('Parent', fig_energy);
  exact_soluitionHandler = plot(ax_energy, ...
    [0, L], ones(1, 2)*get_energy_exact(kbT, J), ...
    'LineWidth', 1.5);
  hold(ax_energy, 'on');
  
  %ylim(ax_energy, [-2, 0]);
  xlim(ax_energy, [0, L]);
  xlabel(ax_energy, 'MCS');
  ylabel(ax_energy, 'E / N');
  
 
  
for l = 1:L
  r0 = ceil(rand(N^2,2)*N);    % site (i, j) on NxN lattice
  rn = mod(r0 - 2,N) + 1;      % for coordinates (i-1, j-1) ; and periodic
  rp = mod(r0,N) + 1;          % for coordinates (i+1, j+1) ; and periodic
  r = rand(N^2,1);
  
  % perform single MCS (Monte Carlo Sweep)
  % we will move randomly from one site to the next site
  % r0(n,:) is random site location.
  for n = 1:N^2
    % (i-1, j), (i+1, j), (i, j-1), (i, j+1)
    % the possible values of the sum nearest neighbours <i, j>
    % of s_i s_j are [-4, -2, 0, 2, 4]
    % division by 2 and addition by 3 are index shifters.
    % Instead of multiply J*s_i*s_j by 2 and divided by 4,
    % we simplified it by dividing by 2, See (Newman and Barkema, 1999) Eq.
    % (3.10)
    i_state = ((x_config_spin(rn(n,1), r0(n,2)) ...
        + x_config_spin(rp(n,1), r0(n,2)) ...
        + x_config_spin(r0(n,1), rn(n,2)) ...
        + x_config_spin(r0(n,1), rp(n,2))) ...
      * x_config_spin(r0(n,1),r0(n,2))/2 + 3);
    
    if r(n) < p(i_state)
      x_config_spin(r0(n,1),r0(n,2)) = -x_config_spin(r0(n,1),r0(n,2)); % flip the spin
    end
  end
  
  imagesc(ax_ising, x_config_spin);
  axis(ax_ising, 'equal', 'off');
  title(ax_ising, sprintf('k_BT=%.1f; samples=(%d/%d)', kbT, l, L));
  
  % -- update data in ax_magnetization
  if mod(l-1, 100) == 0
    M = sum(x_config_spin(:))/N^2;
    scatter(ax_magnetization, l, abs(M), 'ro', 'LineWidth', 1.5);
    title(ax_magnetization, sprintf('k_BT=%.1f, M=%.2f', kbT, abs(M)));
  end
  
  % -- update data in ax_energy
  %    we only compute the pair (avoiding twice calculation by looking at
  %    the neighbours of (i-1, j) and (i, j-1).
  if mod(l-1, 100) == 0
    e1 = 0;
    for i = 1:N
      for j = 1:N
        j_m1 = mod(j-2, N) + 1;  % j-1
        i_m1 = mod(i-2, N) + 1;  % i-1
        %i_p1 = mod(i, N) + 1;    % i+1
        %j_p1 = mod(j, N) + 1;    % j+1
        e_site = -J * x_config_spin(i, j) * ...
          (x_config_spin(i, j_m1) ...
          + x_config_spin(i_m1, j)); % ...
          %+ x_config_spin(i_p1, j) ...
          %+ x_config_spin(i, j_p1));
        
        e1 = e1 + e_site;
      end
    end
    e1 = e1 / N^2;
    
    scatter(ax_energy, l, e1, 'bo', 'LineWidth', 1.5);
    title(ax_energy, sprintf('k_BT=%.1f, E=%.2f', kbT, e1));
  end
  
  drawnow;
end

fprintf('   Computation is finished, total time: %.2f [s]\n', toc);

%% -- function declarations
function energy_exact = get_energy_exact(kbT, J)
  % See `analytical-solution-of-Ising-model.xopp`
  beta = 1/kbT;
  K = 2*sinh(2*beta*J) ./ (cosh(2*beta*J)^2);
  dK_dbeta = 4*J./cosh(2*beta*J).*(1-2*tanh(2*beta*J).^2);
  
  df_dbeta_fun = @(x) -K*sin(x).^2 ./ ...
    ( sqrt(1-K^2.*sin(x).^2) + (1-K^2.*sin(x).^2) ) * dK_dbeta;
  
  df_dbeta_eval = integral(df_dbeta_fun, 0, pi/2, 'AbsTol', 1e-8) / pi;
  
  energy_exact = -tanh(2*beta*J)*2*J - df_dbeta_eval;
  
end

function magnet_exact = get_magnet_exact(kbT, J)
  if kbT > 2/log(1+sqrt(2))
    magnet_exact = 0;
  else
    magnet_exact = real((1 - sinh(2*J*kbT.^-1).^-4).^(1/8));
  end
end