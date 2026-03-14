% integrate the rate equation 
% constant pulling rate
% plot mean y and f versus the rate

clear all;
%digits(50);

%% --- Load energy barriers data
barriers;

n        = n5;     % number of segments - 1
e0       = e0_5;
y_data   = y5;     % y grid
ebs_data = Es5;    % breaking barrier
ebh_data = Eh5;    % healing barrier
ny = length(y_data);

% Interpolation to define ebs and ebh as continuous functions of y
ebs_func = @(y) interp1(y_data, ebs_data, y, 'linear', 'extrap');
ebh_func = @(y) interp1(y_data, ebh_data, y, 'linear', 'extrap');

% Inverse Langevin approximation (valid for |x|<1)
inv_L = @(x) (3*x - x.^3) ./ (1 - x.^2);

%% --- Base simulation parameters
y_target   = y_data(ny);       % target y
N          = 10000;            % time steps

y = linspace(0, y_target, N+1);  
ft = zeros(1, N+1);     % f = f(y)
dydf = zeros(1, N+1);   % compliance dy/df

for i = 1:N+1
    ebs_y(i) = ebs_func(y(i));  % barrier for breaking at current y
    ebh_y(i) = ebh_func(y(i));  % barrier for healing at current y

    ft(i) = inv_L(y(i)/(n+1));
    dydf(i) = ((ft(i))^(-2) - (sinh(ft(i)))^(-2))*(n+1);
end
dydf(1) = (n+1)/3;      % chain compliance at y = 0

ks = (n+1)*exp(-ebs_y);
kh = exp(-ebh_y);

OM = -40:1:-5;
rates = 10.^OM;
nr = length(rates);

y_mean = zeros(1,nr);
f_mean = zeros(1,nr);

y_mean0 = zeros(1,nr);
f_mean0 = zeros(1,nr);

for j = 1:nr
    rate = rates(j);
    total_time = y_target / rate;  % total time
    dt         = total_time / N;   % time step
    dy         = rate * dt;

    P = ones(1, N+1);       % probability for the chain to be intact
    Ps = zeros(1, N+1);     % probability for the chain to be broken

    % integrate the rate equation by backward Euler
    for i = 2:N+1
        P(i) = (P(i-1)+kh(i)*dt)/(1+(ks(i)+kh(i))*dt);
        Ps(i) = (Ps(i-1)+ks(i)*dt)/(1+(ks(i)+kh(i))*dt);
    end

    rhos = P.*ks/rate;      % probability density for scission
    rhoh = Ps.*kh/rate;     % probability density for healing
    rho = rhos-rhoh;

    rhos_f = rhos.*dydf;    % probability density for scission in terms of force
    rhoh_f = rhoh.*dydf;    % probability density for healing in terms of force
    rho_f = rhos_f-rhoh_f;

    % integrate the rate constant for healing from t(N-i) to t(N+1)
    kh_int = zeros(1, N+1);
    for i = 0:N-1
        kh_int(N-i) = kh_int(N-i+1) + (kh(N-i)+kh(N-i+1))/2*dt;
    end

    P_nh = exp(-kh_int);            % probability for no healing after t
    rhos_last = rhos.*P_nh;         % probability density for the last scission event
    rhos_f_last = rhos_f .* P_nh;   % probability density for the last scission event in term of force

    % calculate the average y and f for the last scission events
    for i = 1:N
        y_mean(j) = y_mean(j) + (y(i)*rhos_last(i)+y(i+1)*rhos_last(i+1))/2*dy;
        f_mean(j) = f_mean(j) + (ft(i)*rhos_f_last(i)+ft(i+1)*rhos_f_last(i+1))/2*(ft(i+1)-ft(i));

        y_mean0(j) = y_mean0(j) + (y(i)*rho(i)+y(i+1)*rho(i+1))/2*dy;
        f_mean0(j) = f_mean0(j) + (ft(i)*rho_f(i)+ft(i+1)*rho_f(i+1))/2*(ft(i+1)-ft(i));
    end
end

figure;
semilogx(rates, y_mean/n, 'o-');
%hold on;
%semilogx(rates, y_mean0/n, '^--');
figure;
loglog(rates, f_mean, 'o-');
%hold on;
%loglog(rates, f_mean0, '^--');

fmax = 2.69*e0;
figure;
loglog(rates, f_mean/fmax, 'o-');

