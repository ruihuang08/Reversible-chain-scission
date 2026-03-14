% integrate the rate equation 
% constant pulling rate

clear all;
%digits(50);

%% --- Load energy barriers data
barriers;

y_data   = y3;     % y grid
ebs_data = Es3;    % breaking barrier
ebh_data = Eh3;    % healing barrier
n        = n3;     % number of segments - 1
ny = length(y_data);

% Interpolation to define ebs and ebh as continuous functions of y
ebs_func = @(y) interp1(y_data, ebs_data, y, 'linear', 'extrap');
ebh_func = @(y) interp1(y_data, ebh_data, y, 'linear', 'extrap');

% Inverse Langevin approximation (valid for |x|<1)
inv_L = @(x) (3*x - x.^3) ./ (1 - x.^2);

%% --- Base simulation parameters
rate       = 1e-9;            % y(t) = rate * t
%y_target   = n;              % target y
y_target   = y_data(ny);       % target y
total_time = y_target / rate;  % total time
N          = 10000;            % time steps
dt         = total_time / N;   % time step

y = rate * linspace(0, total_time, N+1);   % y(t) = rate * t
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

P = ones(1, N+1);       % probability for the chain to be intact
Ps = zeros(1, N+1);     % probability for the chain to be broken
P_noh = ones(1, N+1);   % probability for the chain to be intact when healing is off

% integrate the rate equation by backward Euler
for i = 2:N+1
    P(i) = (P(i-1)+kh(i)*dt)/(1+(ks(i)+kh(i))*dt);
    Ps(i) = (Ps(i-1)+ks(i)*dt)/(1+(ks(i)+kh(i))*dt);
    P_noh(i) = P_noh(i-1)/(1+ks(i)*dt); % no healing
end

rhos = P.*ks/rate;      % probability density for scission
rhoh = Ps.*kh/rate;     % probability density for healing
rhos_noh = P_noh.*ks/rate;      % probability density for scission

rhos_f = rhos.*dydf;    % probability density for scission in terms of force
rhoh_f = rhoh.*dydf;    % probability density for healing in terms of force
rhos_f_noh = rhos_noh.*dydf;    % probability density for scission in terms of force

%figure;
%plot(y, P);
%hold on;
%plot(y, P_noh);
%figure;
%plot(y, rhos, 's');
%hold on;
%plot(y, rhoh, '^');
%plot(y, rhos-rhoh);

%rhos_noh = P_noh.*ks/rate;
%plot(y, rhos_noh, 'o');

% integrate the probability density for healing and scission
Ph_tot = 0;
Ps_tot = 0;
for i = 1:N
    Ph_tot = Ph_tot + (rhoh(i)+rhoh(i+1))/2*rate*dt;
    Ps_tot = Ps_tot + (rhos(i)+rhos(i+1))/2*rate*dt;
end

% cumulative probability for healing at y > yi
Phi = zeros(1,N+1);
Phi(1) = Ph_tot;
for i = 2:N+1
    Phi(i) = Phi(i-1) - (rhoh(i)+rhoh(i-1))/2*rate*dt;
end

% integrate the rate constant for healing from t(N-i) to t(N+1)
kh_int = zeros(1, N+1);
for i = 0:N-1
    kh_int(N-i) = kh_int(N-i+1) + (kh(N-i)+kh(N-i+1))/2*dt;
end

P_nh = exp(-kh_int);            % probability for no healing after t
rhos_last = rhos.*P_nh;         % probability density for the last scission event
rhos_f_last = rhos_f .* P_nh;   % probability density for the last scission event in term of force

%calculate the area under the curve rhos_last
% both A and Af should be close to 1
A = 0;
Af = 0;
for i = 1:N
    A = A + (rhos_last(i)+rhos_last(i+1))/2*dt*rate;
    Af = Af + (rhos_f_last(i)+rhos_f_last(i+1))/2*(ft(i+1)-ft(i));
end

figure;
plot(y, rhos_last, 'LineWidth',2);
hold on;
plot(y, rhos, ':', 'LineWidth',2);
plot(y, rhoh, '-.', 'LineWidth',2);
plot(y, rhos-rhoh, '--', 'LineWidth',2);
plot(y, rhos_noh, ':', 'LineWidth',2);

figure;
plot(ft, rhos_f_last, 'LineWidth',2);
hold on;
plot(ft, rhos_f, ':', 'LineWidth',2);
plot(ft, rhoh_f, '-.', 'LineWidth',2);
plot(ft, rhos_f-rhoh_f, '--', 'LineWidth',2);
plot(ft, rhos_f_noh, ':', 'LineWidth',2);

% calculate the average y and f for the last scission events
y_mean = 0;
f_mean = 0;

for i = 1:N
    y_mean = y_mean + (y(i)*rhos_last(i)+y(i+1)*rhos_last(i+1))/2*(y(i+1)-y(i));
    f_mean = f_mean + (ft(i)*rhos_f_last(i)+ft(i+1)*rhos_f_last(i+1))/2*(ft(i+1)-ft(i));
end


