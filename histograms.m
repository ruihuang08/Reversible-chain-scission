% run KMC simulations at a constant pulling rate 
% create histograms for the rupture forces and end-to-end distance at the
% final scission event
% compare to the prediction by the rate equation

clear all;

%% --- Load data from the Excel file
%T = readtable('barriers.xlsx');
%y_data   = T.y_n100_e50;      % y grid
%ebs_data = T.ebs_n100_e50;    % breaking barrier
%ebh_data = T.ebh_n100_e50;    % healing barrier

barriers;
y_data   = y3;     % y grid
ebs_data = Es3;    % breaking barrier
ebh_data = Eh3;    % healing barrier

% Interpolation to define ebs and ebh as continuous functions of y
ebs_func = @(y) interp1(y_data, ebs_data, y, 'linear', 'extrap');
ebh_func = @(y) interp1(y_data, ebh_data, y, 'linear', 'extrap');

% Inverse Langevin approximation (valid for |x|<1)
inv_L = @(x) (3*x - x.^3) ./ (1 - x.^2);

%% --- Base simulation parameters
n          = 100;              % n+1 is the number of links in the chain
y_target   = n;                % target y
rate       = 1e-19;            % y(t) = rate * t
total_time = y_target / rate;  % total time
N          = 10000;            % time steps
dt         = total_time / N;   % time step
healing_switch = 'on';         % allow healing before final scission

%% integrate the rate equation
yt = rate * linspace(0, total_time, N+1);   % y(t) = rate * t
ft = zeros(1, N+1);     % f = f(y)
dydf = zeros(1, N+1);   % compliance dy/df

for i = 1:N+1
    ebs_y(i) = ebs_func(yt(i));  % barrier for breaking at current y
    ebh_y(i) = ebh_func(yt(i));  % barrier for healing at current y

    ft(i) = inv_L(yt(i)/(n+1));
    dydf(i) = ((ft(i))^(-2) - (sinh(ft(i)))^(-2))*(n+1);
end
dydf(1) = (n+1)/3;      % chain compliance at y = 0

ks = (n+1)*exp(-ebs_y); % rate constance for scission
kh = exp(-ebh_y);       % rate constant for healing
P = ones(1, N+1);       % probability for the chain to be intact at t
Ps = zeros(1, N+1);     % probability for the chain to be broken at t

% integrate the rate equation by backward Euler
for i = 2:N+1
    P(i) = (P(i-1)+kh(i)*dt)/(1+(ks(i)+kh(i))*dt);
    Ps(i) = (Ps(i-1)+ks(i)*dt)/(1+(ks(i)+kh(i))*dt);
end

rhos = P.*ks/rate;      % probability density for scission
rhoh = Ps.*kh/rate;     % probability density for healing

rhos_f = rhos.*dydf;    % probability density for scission in terms of force
rhoh_f = rhoh.*dydf;    % probability density for healing in terms of force

% integrate the rate constant for healing from t(N-i) to t(N+1)
kh_int = zeros(1, N+1);
for i = 0:N-1
    kh_int(N-i) = kh_int(N-i+1) + (kh(N-i)+kh(N-i+1))/2*dt;
end

P_nh = exp(-kh_int);            % probability for not healing between t(i) to t(N+1)
rhos_last = rhos .* P_nh;       % probability density for the last scission event
rhos_f_last = rhos_f .* P_nh;   % probability density for the last scission event in term of force

%calculate the area under the pdf curve rhos_last
A = 0;
Af = 0;
for i = 1:N
    A = A + (rhos_last(i)+rhos_last(i+1))/2*dt*rate;
    Af = Af + (rhos_f_last(i)+rhos_f_last(i+1))/2*(ft(i+1)-ft(i));
end

rhos_last_normalized = rhos_last/A;         % this is unnecessary, becasue A = 1
rhos_f_last_normalized = rhos_f_last/Af; 

%% --- Multi-run settings
M = 500;                       % number of trajectories
scission_force = NaN(M,1);
scission_y     = NaN(M,1);

h = waitbar(0,'Running simulations...');  % progress bar

for m = 1:M
    [y, bond_state, force] = run_single_trajectory(N, dt, rate, n, inv_L, ebs_func, ebh_func, healing_switch);

    % Find all 1->0 transitions
    trans_idx = find(bond_state(1:end-1)==1 & bond_state(2:end)==0);

    if ~isempty(trans_idx)
        idx = trans_idx(end); % final scission index
        scission_force(m) = force(idx);
        scission_y(m)     = y(idx);
    end

    % update progress bar
    if mod(m,10)==0 || m==M
        waitbar(m/M,h,sprintf('Simulation %d of %d',m,M));
    end
end
close(h);

%theory = readtable('rho_data.xlsx','Sheet','all_data','VariableNamingRule','preserve');

% Extract columns for rate = 1e-16 (use exact Excel headers)
%y_col    = theory.("y_1e-16");
%f_col    = theory.("f_1e-16");
%rho_y_th = theory.("rho_y_1e-16");
%rho_f_th = theory.("rho_f_1e-16");

% Clean NaNs (from padding in the Excel export)
%valid_y = ~isnan(y_col) & ~isnan(rho_y_th);
%valid_f = ~isnan(f_col) & ~isnan(rho_f_th);

%y_col    = y_col(valid_y);
%rho_y_th = rho_y_th(valid_y);
%f_col    = f_col(valid_f);
%rho_f_th = rho_f_th(valid_f);



%% --- Histogram of scission forces with overlay
figure;
h1 = histogram(scission_force(~isnan(scission_force)), ...
    'BinMethod','fd','LineWidth',1.5);
hold on;

% Scale rho_f by counts: M * binWidth
binWidthF = h1.BinWidth;
%rho_f_scaled = rhos_f_last_normalized * M * binWidthF;
rho_f_scaled = rhos_f_last * M * binWidthF;

rho_f_scaled2 = (rhos_f-rhoh_f) * M * binWidthF;

plot(ft, rho_f_scaled, 'r-', 'LineWidth', 2);
plot(ft, rho_f_scaled2, 'r--', 'LineWidth', 2);

xlabel('$\bar{f}_p$','Interpreter','latex','FontSize',24);
ylabel('Count','Interpreter','latex','FontSize',24);
set(gca,'FontSize',20); grid on;
legend('KMC','\rho_{s,last}', '\rho_s-\rho_h');

%% --- Histogram of scission y with overlay
figure;
h2 = histogram(scission_y(~isnan(scission_y)), ...
    'BinMethod','fd','LineWidth',1.5);
hold on;

% Scale rho_y by counts: M * binWidth
binWidthY = h2.BinWidth;
%rho_y_scaled = rhos_last_normalized * M * binWidthY;
rho_y_scaled = rhos_last * M * binWidthY;

rho_y_scaled2 = (rhos-rhoh) * M * binWidthY;

plot(yt, rho_y_scaled, 'r-', 'LineWidth', 2);
plot(yt, rho_y_scaled2, 'r--', 'LineWidth', 2);

xlabel('$\bar{y}_p$','Interpreter','latex','FontSize',24);
ylabel('Count','Interpreter','latex','FontSize',24);
set(gca,'FontSize',20); grid on;
legend('KMC','\rho_{s,last}', '\rho_s-\rho_h');

%% ===== Local function =====
function [y, bond_state, force] = run_single_trajectory(N, dt, rate, n, inv_L, ebs_func, ebh_func, healing_switch)
    % Preallocate
    y = rate * linspace(0, dt*(N-1), N);   % y(t) = rate * t
    bond_state = zeros(1, N);
    bond_state(1) = 1;                     % start healed
    force = zeros(1, N);

    for i = 1:N-1
        ebs_y = ebs_func(y(i));  % barrier for breaking at current y
        ebh_y = ebh_func(y(i));  % barrier for healing at current y

        switch healing_switch
            case 'on'
                if bond_state(i) == 1
                    p = (n+1) * dt * exp(-ebs_y);
                    bond_state(i+1) = rand() >= p;
                else
                    p = dt * exp(-ebh_y);
                    bond_state(i+1) = rand() < p;
                end
            case 'off'
                if bond_state(i) == 1
                    p = (n+1) * dt * exp(-ebs_y);
                    bond_state(i+1) = rand() >= p;
                else
                    bond_state(i+1) = 0; % stays broken
                end
            otherwise
                error('Invalid healing_switch. Use ''on'' or ''off''.');
        end

        % Compute force at step i using the *state at step i*
        if bond_state(i) == 1
            x = y(i) / (n+1);
            if abs(x) < 1
                force(i) = inv_L(x);
            else
                force(i) = NaN;
            end
        else
            force(i) = 0;
        end
    end

    % Force at last index
    if bond_state(N) == 1
        xN = y(N) / (n+1);
        if abs(xN) < 1
            force(N) = inv_L(xN);
        else
            force(N) = NaN;
        end
    else
        force(N) = 0;
    end
end
