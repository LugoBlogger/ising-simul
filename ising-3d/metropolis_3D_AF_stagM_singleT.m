% metropolis_3D_AF_stagM_singleT.m
%
%  This script simualtes antiferromagnetic Ising model in a cubic lattice
%  We do not have an exact solution for the magnetization and we only test
%  the equilibration.
%
%
%
% To do list
% [x] Improve scatter3 plot to update only the CData
%     without redraw from the beginning

clear;

%% -- user-defined input
N = 10;
L = 5000;   % number of samples after one MCS (Monte Carlo Sweep)
J = -1;
kb = 1;
T_i = 4.2;   % Tc = ??? (in 3D cubic lattice)


kbT = kb*T_i;   % we set k_B to 1

% initial configuration of spins are random
x_config_spin = ceil(rand(N, N, N)*2)*2 - 3;   % (0, 2) -> {1, 2} -> {2, 4} -> {-1, 1};

% probability that guides each flip is allowed or not
% this probability is determined from acceptance ratio
p = [exp(12*J/kbT), ...   % Delta H = -12 abs(J) [allowed with probability]
     exp(8*J/kbT), ...    % Delta H = -8 abs(J) [allowed with probability]
     exp(4*J/kbT), ...    % Delta H = -4 abs(J) [allowed with probability]
     1, ...               % Delta H = 0 -> p(0) = exp(0) = 1
     1, ...               % Delta H = 4 abs(J)
     1, ...               % Delta H = 8 abs(J)
     1];                  % Delta H = 12 abs(J)
   
yMultiplier = 1;
zMultiplier = 8;
markerSize = 20;
[xMesh, yMesh, zMesh] = meshgrid(1:N, 1:yMultiplier:yMultiplier*N, 1:zMultiplier:zMultiplier*N);
xData = xMesh(:);
yData = yMesh(:);
zData = zMesh(:);

%% -- run Monte Carlo simulation for 2D Ising model
tic;

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
  ylim(ax_magnetization, [-0.1, 1.1]);
  xlim(ax_magnetization, [0, L]);
  ylabel(ax_magnetization, 'stag. mag. / N');
  
for l = 1:L
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
  
  set(scatterHandler, 'CData', get_spinColorData(x_config_spin));
  title(ax_ising, sprintf('k_BT=%.1f; samples=(%d/%d)', kbT, l, L));
  
  
  % -- update data in ax_magnetization
  if mod(l-1, 100) == 0
    M = get_spinSum_staggered(x_config_spin, N)/N^3;
    scatter(ax_magnetization, l, abs(M), 'ro', 'LineWidth', 1.5);
    title(ax_magnetization, sprintf('k_BT=%.1f, M=%.2f', kbT, abs(M)));
  end
  
  
  drawnow;
end
  
   
%% -- function declaration
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

function spinColorData = get_spinColorData(x_config_spin)
%   colorSpace = [1.0, 0.0, 0.0; 
%                 0.0, 0.0, 1.0];
  colorSpace = [0.2422, 0.1504, 0.6603;    % parula (https://stackoverflow.com/a/60007513)
              0.9769, 0.9839, 0.0805];

  idx_x_config_spin = (x_config_spin+3)/2;
  idx_x_config_spin = int8(idx_x_config_spin(:));
  spinColorData = colorSpace(idx_x_config_spin,:);
end