% metropolis_3D_AF_all_T.m
%
%  This script is based on metropolis_AF_all_T.m. To calculate the
%  magnetization, we used a staggered magnetization
%
% Logs
%
%
%
%
clear;

%% -- user-defined input
N = 10;
J = -1;     % J=-1 (anti-ferromagnetic)
kb = 1;
Tmin = 3;     % (Sosin et al., 2015) - Computational analysis of 3D Ising model using metropolis
Tmax = 6;

N_T = 50;       % number of grid points along temperature axis; default 25

L_avg = 4000;   % number of lastest states that we used in average computation
              % default 100

L_toHighT = 4000;   % number of states that we have on each temperature setting
                  % default 1000
L_toLowT = 5000;    % for the direction from high temperature to low temperature
                  % default 8000
isToHighT = false; 

T_arr = linspace(Tmin, Tmax, N_T);
[kbT, dT, x_config_spin, init_energy, L] = get_init_setting(N, kb, J, T_arr, ...
  L_toHighT, L_toLowT, isToHighT);


Q = L * N_T;    % total number of generated samples for all temperatures
spinSum_arr = zeros(L, 1);  % array to store sum of spin
energy_arr = zeros(L, 1);  % this will store the sampled energies for each temperature

%% -- compute and simulate expectation value of magnetization
computeObservables = [1, 0, 1];     % [magnetization, susceptibility, specific heat]
simulate_IsingModel(J, Tmin, Tmax, kb, kbT, N, L, L_avg, ...
  x_config_spin, init_energy, dT, Q, spinSum_arr, energy_arr, computeObservables)


%% -- function declaration
function simulate_IsingModel(J, Tmin, Tmax, kb, kbT, N, L, L_avg, ...
  x_config_spin, energy_each_state, dT, Q, spinSum_arr, energy_arr, computeObservables)

  tic;
  
  %beta_arr = 1./(kb*Tf);
  yMultiplier = 1;
  zMultiplier = 8;
  markerSize = 20;
  [xMesh, yMesh, zMesh] = meshgrid(1:N, 1:yMultiplier:yMultiplier*N, 1:zMultiplier:zMultiplier*N);
  xData = xMesh(:);
  yData = yMesh(:);
  zData = zMesh(:);
  
  
  % -- ising simulation figure handle
  fig_ising = figure('Units', 'centimeters', 'Position', [3, 1.5, 8, 16], 'Color', 'w');
    ax_ising = axes('Parent', fig_ising, 'Position', [0.02, 0.02, 0.95, 0.92]);
    
    scatterHandler = scatter3(xData, yData, zData, markerSize, ...
    get_spinColorData(x_config_spin), 'filled', 'MarkerFaceAlpha', 0.7, ...
    'MarkerEdgeColor', 'none');
    axis(ax_ising, 'equal', 'off');
    view(ax_ising, [60, 30]);
    
  % -- magnetization curve figure handle
  fig_magnetization = figure('Units', 'centimeters', 'Position', [16, 3, 20, 12], 'Color', 'w');
    ax_magnetization = axes('Parent', fig_magnetization);
    hold(ax_magnetization, 'on');
    xlim(ax_magnetization, [Tmin, Tmax]);
    ylim(ax_magnetization, [-0.1, 1.1]);
    
  % -- specific heat curve figure handle
  fig_spec_heat = figure('Units', 'centimeters', 'Position', [17, 4, 20, 12], 'Color', 'w');
    ax_spec_heat = axes('Parent', fig_spec_heat);
    
    xlim(ax_spec_heat, [Tmin, Tmax]);
    %ylim(ax_spec_heat, [-0.1, 2.1]);
    hold(ax_spec_heat, 'on');
    
  % -- susceptibility curve figure handle
  fig_suscep = figure('Units', 'centimeters', 'Position', [18, 5, 20, 12], 'Color', 'w');
    ax_suscep = axes('Parent', fig_suscep);
    
    xlim(ax_suscep, [Tmin, Tmax]);
    hold(ax_suscep, 'on');  
  
    
    
  for q = 1:Q
    if mod(q-1, L) == 0
      message = sprintf('   kbT: %.3f', kbT);
      fprintf(message);    
    end

    % probability that guides each flip is allowed or not
    % this probability is determined from acceptance ratio
    p = [exp(12*J/kbT), ...   % Delta H = -12 abs(J) [allowed with probability]
         exp(8*J/kbT), ...    % Delta H = -8 abs(J) [allowed with probability]
         exp(4*J/kbT), ...    % Delta H = -4 abs(J) [allowed with probability]
         1, ...               % Delta H = 0 -> p(0) = exp(0) = 1
         1, ...               % Delta H = 4 abs(J)
         1, ...               % Delta H = 8 abs(J)
         1];                  % Delta H = 12 abs(J)
    
    r0 = ceil(rand(N^3, 3)*N);        % site (i, j, k) on N x N x N latticc
    rn = mod(r0 - 2, N) + 1;          % for coordinates (i-1, j-1, k-1) ; and periodic
    rp = mod(r0, N) + 1;              % for coordinates (i+1, j+1, k+1) ; and periodic
    r = rand(N^3, 1);
    
    % perform single MCS (Monte Carlo Sweep)
    % we will move randomly from one site to the next site
    % r0(n, :) is random site location.
    for n = 1:N^3
      % (i-1, j, k), (i+1, j, k), 
      % (i, j-1, k), (i, j+1, k),
      % (i, j, k-1), (i, j, k+1)
      % the possible values of the sum nearest neighbours <K, L>
      % of s_K s_L are [-6, -4, -2, 0, 2, 4, 6]
      % division by 2 and addition by 4 are index shifters.
      % Instead of multiply J*s_K*s_L (K and L are coordinates in a cubic
      % lattice) by 2 and divided by 4, we simplified it by dividing by 2,
      % See (Newman and Barkema, 1999) Eq. (3.10) for a square lattice
      i_state = ((x_config_spin(rn(n, 1), r0(n, 2), r0(n, 3)) ...
          + x_config_spin(rp(n, 1), r0(n, 2), r0(n, 3)) ...
          + x_config_spin(r0(n, 1), rn(n, 2), r0(n, 3)) ...
          + x_config_spin(r0(n, 1), rp(n, 2), r0(n, 3)) ...
          + x_config_spin(r0(n, 1), r0(n, 2), rn(n, 3)) ...
          + x_config_spin(r0(n, 1), r0(n, 2), rp(n, 3))) ...
        * x_config_spin(r0(n, 1), r0(n, 2), r0(n, 3))/2 + 4);

      if r(n) < p(i_state)
        x_config_spin(r0(n, 1), r0(n, 2), r0(n, 3)) ...
          = -x_config_spin(r0(n, 1), r0(n, 2), r0(n, 3));  % flip the spin 
      end
    end
    
    
    % -- update ising plot
    set(scatterHandler, 'CData', get_spinColorData(x_config_spin));
    title(ax_ising, sprintf('k_BT=%.1f;  samples=(%d/%d)', kbT, mod(q, L)+1, L));
  
  
    % -- update energy_arr
    energy_each_state = get_energy(x_config_spin, J, N);
    energy_arr = [energy_arr(2:end); energy_each_state];   % First In First Out
  
    
    % -- update magnetization_arr (staggered magnetization)
    spinSum_each_state = get_spinSum_staggered(x_config_spin, N)/(N^3);
    spinSum_arr = [spinSum_arr(2:end);  spinSum_each_state];  % First In First Out
  
  
    % -- update data in ax_magnetization
    if mod(q, L) == 0
      M = mean(spinSum_arr(L-(L_avg-1):L)); 
      scatter(ax_magnetization, kbT, abs(M), 'ro', 'LineWidth', 1.5);
      title(ax_magnetization, sprintf('k_BT=%.1f, M=%.2f', kbT, abs(M)));
    end
    
    % -- update data in ax_spec_heat
    if mod(q, L) == 0
      C = kb * var(energy_arr(L-(L_avg-1):L))/((N*kbT)^2); % See (Newman and Barkema, 1999) Eq. (1.37)
      
      scatter(ax_spec_heat, kbT, C, 'bo', 'LineWidth', 1.5);
      title(ax_spec_heat, sprintf('k_BT=%.1f, C=%.2f', kbT, C));
    end
    
    % -- update data in ax_suscep
    if mod(q, L) == 0
      chi_suscep = N*var(spinSum_arr(L-(L_avg-1):L))/kbT;
      
      scatter(ax_suscep, kbT, chi_suscep, 'mo', 'LineWidth', 1.5);
      title(ax_suscep, sprintf('k_BT=%.1f, chi=%.2f', kbT, chi_suscep));
    end

    drawnow;
  
    
    if mod(q, L) == 0
      fprintf(repmat('\b', 1, length(message)));

      kbT = kbT + dT;
    end
  end
  
  fprintf('   Computation is finished, total time: %.2f [s]\n', toc);
end

function spinSum_each_state = get_spinSum_staggered(x_config_spin, N)
  spinSum_each_state = 0;
  checkboardMatrix = ones(N, N, N);
  checkboardMatrix(1:2:N, 2:2:N, 1:2:N) = -1;
  checkboardMatrix(2:2:N, 1:2:N, 1:2:N) = -1;
  
  checkboardMatrix(1:2:N, 1:2:N, 2:2:N) = -1;
  checkboardMatrix(2:2:N, 2:2:N, 2:2:N) = -1;
 
  for i = 1:N
    for j = 1:N
      for k = 1:N
        spinSum_each_state = spinSum_each_state + checkboardMatrix(i, j, k)*x_config_spin(i, j, k);
      end
    end
  end
  
end

function energy_each_state = get_energy(x_config_spin, J, N)
  % we only compute the pair (avoiding twice calculation by looking at
  % the neighbours of (i-1, j) and (i, j-1).
  % the following procedures are similar to the calculation of init_energy
  energy_each_state = 0;
  for i = 1:N
    for j = 1:N
      for k = 1:N
        k_m1 = mod(k-2, N) + 1;  % k-1
        j_m1 = mod(j-2, N) + 1;  % j-1
        i_m1 = mod(i-2, N) + 1;  % i-1
        e_site = -J * x_config_spin(i, j, k) * ...
          ( x_config_spin(i, j, k_m1) ...
          + x_config_spin(i, j_m1, k) ...
          + x_config_spin(i_m1, j, k)); 

        energy_each_state = energy_each_state + e_site;
      end
    end
  end
end

function [kbT, dT, x_config_spin, init_energy, L] = get_init_setting(N, kb, J, T_arr, ...
  L_toHightT, L_toLowT, isToHighT)
  if isToHighT 
    kbT = kb*T_arr(1);
    dT = T_arr(2) - T_arr(1);
    x_config_spin = get_checkerboardInCubic(N);
    L = L_toHightT;
  else
    kbT = kb*T_arr(end);
    dT = T_arr(1) - T_arr(2);
    x_config_spin = ceil(rand(N, N, N)*2)*2 - 3;   % (0, 2) -> {1, 2} -> {2, 4} -> {-1, 1};
    L = L_toLowT;
  end
  
  % get init_energy
  % we only compute the pair (avoiding twice calculation by looking at
  % the neighbours of (i-1, j) and (i, j-1).
  init_energy = 0;
  for i = 1:N
    for j = 1:N
      for k = 1:N
        k_m1 = mod(k-2, N) + 1;  % k-1
        j_m1 = mod(j-2, N) + 1;  % j-1
        i_m1 = mod(i-2, N) + 1;  % i-1
        e_site = -J * x_config_spin(i, j, k) * ...
          ( x_config_spin(i   , j   , k_m1) ...
          + x_config_spin(i   , j_m1, k) ...
          + x_config_spin(i_m1, j   , k)); 

        init_energy = init_energy + e_site;
      end
    end
  end
  
end

function spinColorData = get_spinColorData(x_config_spin)
%   colorSpace = [1.0, 0.0, 0.0; 
%                 0.0, 0.0, 1.0];
  colorSpace = [0.2422, 0.1504, 0.6603;    % parula (https://stackoverflow.com/a/60007513)
              0.9769, 0.9839, 0.0805];

  idx_x_config_spin = (x_config_spin+3)/2;
  idx_x_config_spin = int8(idx_x_config_spin(:));
  spinColorData = colorSpace(idx_x_config_spin,:);
end

function checkboardInCubic = get_checkerboardInCubic(N)
  checkboardInCubic = ones(N, N, N);
  checkboardInCubic(1:2:N, 2:2:N, 1:2:N) = -1;
  checkboardInCubic(2:2:N, 1:2:N, 1:2:N) = -1;
  
  checkboardInCubic(1:2:N, 1:2:N, 2:2:N) = -1;
  checkboardInCubic(2:2:N, 2:2:N, 2:2:N) = -1;
end
  