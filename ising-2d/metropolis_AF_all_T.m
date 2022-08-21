% metropolis_AF_all_T.m
%
%  This script is based on metropolis_2_v2.m. To calculate the
%  magnetization, we used a staggered magnetization
%
% Logs
%
%
%
% To do:
% - Revise the forumla update of magnetization and energy for each
%   flipping. We do this if the number of sites is larger than the number
%   of generated samples.

clear;

%% -- user-defined input
N = 40;     % default 50
J = -1;      % J=-1 (anti-ferromagnetic)
kb = 1;
Tmin = 0; 
Tmax = 5; 
%Tmax = 2.551;


N_T = 50;     % number of grid point along temperature axis; default 25
%N_T = 26;

L_avg = 4000;   % number of lastest states that we used in average computation
              % default 100

L_toHighT = 4000;   % number of states that we have on each temperature setting
                  % default 1000
L_toLowT = 5000;    % for the direction from high temperature to low temperature
                  % default 8000
isToHighT = false;              

Tc = 2*abs(J)/(kb*log(1+sqrt(2)));    % exact solution
Tf = 0:0.001:Tc;         % the range of temperature to calculate exact solution of magnetization

T_arr = linspace(Tmin, Tmax, N_T);
[kbT, dT, x_config_spin, init_energy, L] = get_init_setting(N, kb, J, T_arr, ...
  L_toHighT, L_toLowT, isToHighT);

Q = L * N_T;    % total number of generated samples for all temperatures
spinSum_arr = zeros(L, 1);  % array to store sum of spin
energy_arr = zeros(L, 1);  % this will store the sampled energies for each temperature

%% -- compute and simulate expectation value of magnetization
computeObservables = [1, 0, 1];     % [magnetization, susceptibility, specific heat]
simulate_IsingModel(Tf, J, Tc, Tmin, Tmax, kb, kbT, N, L, L_avg, ...
  x_config_spin, init_energy, dT, Q, spinSum_arr, energy_arr, computeObservables)


%% -- function declarations
function simulate_IsingModel(Tf, J, Tc, Tmin, Tmax, kb, kbT, N, L, L_avg, ...
  x_config_spin, energy_each_state, dT, Q, spinSum_arr, energy_arr, computeObservables)

  tic;
  
  beta_arr = 1./(kb*Tf);
  %spinSum_each_state = get_spinSum_staggered(x_config_spin, N);

  % -- ising simulation figure handle
  fig_ising = figure('Units', 'centimeters', 'Position', [3, 3, 12, 12], 'Color', 'w');
    ax_ising = axes('Parent', fig_ising);

  % -- magnetization curve figure handle
  fig_magnetization = figure('Units', 'centimeters', 'Position', [16, 3, 20, 12], 'Color', 'w');
    ax_magnetization = axes('Parent', fig_magnetization);
    exact_solutionHandler = plot(ax_magnetization, Tf, get_magnet_exact(beta_arr, J), ...
      'LineWidth', 1.5);


    hold(ax_magnetization, 'on');
      plot(ax_magnetization, [Tc, Tmax], [0, 0], 'LineWidth', 1.5, ...
        'Color', get(exact_solutionHandler, 'Color'));
      %magnetization_plot_handler = scatter(ax_magnetization, [], [], 'ro', ...
      %  'LineWidth', 1.5);
    %hold(ax_magnetization, 'off');

    xlim(ax_magnetization, [Tmin, Tmax]);
    ylim(ax_magnetization, [-0.1, 1.1]);
  
  % -- specific heat curve figure handle
  fig_spec_heat = figure('Units', 'centimeters', 'Position', [17, 4, 20, 12], 'Color', 'w');
    ax_spec_heat = axes('Parent', fig_spec_heat);
    
    xlim(ax_spec_heat, [Tmin, Tmax]);
    ylim(ax_spec_heat, [-0.1, 2.1]);
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

    % probability that guides each flip is allowed or not (different from
    % ferromagnetic ising model, configuration idx 0 and 1 have allowed
    % transition probability)
    p = [exp(8*J/kbT), ...   % Delta H = 8J [allowed with probability]
         exp(4*J/kbT), ...   % Delta H = 4J [allowed with probability]
         1, ...              % Delta H = 0 -> p(0) = exp(0) = 1
         1, ...  % Delta H = -4J 
         1];     % Delta H = -8J 

    r0 = ceil(rand(N^2, 2) * N);   % site (i, j) on NxN lattice
    rn = mod(r0 - 2, N) + 1;       % for coordinates (i-1, j-1)
    rp = mod(r0, N) + 1;           % for coordinates (i+1, j+1)
    r = rand(N^2,1);    

    % we will move randomly from one site to the next site
    % r0(n,:) is random site location.
    for n = 1:N^2
      % (i-1, j), (i+1, j), (i, j-1), (i, j+1)
      % the possible values of the sum nearest neighbours <i, j>
      % of s_i s_j are [-4, -2, 0, 2, 4]
      % division by 2 and addition by 3 are index shifters.
      i_state = ((x_config_spin(rn(n,1), r0(n,2)) ...
          + x_config_spin(rp(n,1), r0(n,2)) ...
          + x_config_spin(r0(n,1), rn(n,2)) ... 
          + x_config_spin(r0(n,1), rp(n,2))) ...
        * x_config_spin(r0(n,1), r0(n,2))/2 + 3);

      if r(n) < p(i_state)
        x_config_spin(r0(n,1), r0(n,2)) = -x_config_spin(r0(n,1), r0(n,2));
        
        % -- compute the energy of the next sampled state
        %    See (Newman and Barkema, 1999) Eq. 3.11
        %    use the following lines if the number of sites is larger than
        %    the number of iterations
%         deltaE = 2*J*(i_state - 3)*2;
%         energy_each_state = energy_each_state + deltaE;
        
        % -- compute the magnetization of the next sampled state
        %    See (Newman and Barkema, 1999) Eq. 3.14
        %    use the following lines if the number of sites is larger than
        %    the number of iterations
%         deltaM = 2*x_config_spin(r0(n,1), r0(n,2));
%         spinSum_each_state = spinSum_each_state + deltaM;
      end
    end

    % -- update ising plot
    imagesc(ax_ising, x_config_spin, [-1, 1]);
    axis(ax_ising, 'equal', 'off');
    title(ax_ising, sprintf('k_BT=%.1f;  samples=(%d/%d)', kbT, mod(q, L)+1, L));
    
    % -- update energy_arr
    energy_each_state = get_energy(x_config_spin, J, N);
    energy_arr = [energy_arr(2:end); energy_each_state];   % First In First Out
    

    % -- update magnetization_arr (staggered magnetization)
    spinSum_each_state = get_spinSum_staggered(x_config_spin, N)/(N^2);
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
  checkboardMatrix = ones(N, N);
  checkboardMatrix(1:2:N, 2:2:N) = -1;
  checkboardMatrix(2:2:N, 1:2:N) = -1;
 
  for i = 1:N
    for j = 1:N
      spinSum_each_state = spinSum_each_state + checkboardMatrix(i, j)*x_config_spin(i, j);
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
      j_m1 = mod(j-2, N) + 1;  % j-1
      i_m1 = mod(i-2, N) + 1;  % i-1
      e_site = -J * x_config_spin(i, j) * ...
        (x_config_spin(i, j_m1) ...
        + x_config_spin(i_m1, j)); 

      energy_each_state = energy_each_state + e_site;
    end
  end
end

function [kbT, dT, x_config_spin, init_energy, L] = get_init_setting(N, kb, J, T_arr, ...
  L_toHightT, L_toLowT, isToHighT)
  if isToHighT 
    kbT = kb*T_arr(1);
    dT = T_arr(2) - T_arr(1);
    x_config_spin = ones(N, N);
    L = L_toHightT;
  else
    kbT = kb*T_arr(end);
    dT = T_arr(1) - T_arr(2);
    x_config_spin = ceil(rand(N,N)*2)*2 - 3;   % [0, 2] -> {1, 2} -> {2, 4} -> {-1, 1};
    L = L_toLowT;
  end
  
  % get init_energy
  % we only compute the pair (avoiding twice calculation by looking at
  % the neighbours of (i-1, j) and (i, j-1).
  init_energy = 0;
  for i = 1:N
    for j = 1:N
      j_m1 = mod(j-2, N) + 1;  % j-1
      i_m1 = mod(i-2, N) + 1;  % i-1
      e_site = -J * x_config_spin(i, j) * ...
        (x_config_spin(i, j_m1) ...
        + x_config_spin(i_m1, j)); 

      init_energy = init_energy + e_site;
    end
  end
  
end

function magnet_exact = get_magnet_exact(beta, J)
  magnet_exact = real((1 - sinh(2*J*beta).^-4).^(1/8));
end
