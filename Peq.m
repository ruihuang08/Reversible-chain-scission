% equilibrium probability density functions

clear all;

%% --- Load energy barriers data
barriers;

n = n5;
y_data   = y5;     % y grid
ebs_data = Es5;    % breaking barrier
ebh_data = Eh5;    % healing barrier

% Interpolation to define ebs and ebh as continuous functions of y
ebs_func = @(y) interp1(y_data, ebs_data, y, 'linear', 'extrap');
ebh_func = @(y) interp1(y_data, ebh_data, y, 'linear', 'extrap');

% Inverse Langevin approximation (valid for |x|<1)
inv_L = @(x) (3*x - x.^3) ./ (1 - x.^2);

%% --- Base simulation parameters
y_target   = n;              % target y
Ny         = 500;          % time steps
dy         = y_target / Ny;  % y step

y = linspace(0, y_target, Ny+1);  
force = zeros(1, Ny+1);
dydf = zeros(1, Ny+1);   % compliance dy/df

for i = 1:Ny+1
    ebs_y(i) = ebs_func(y(i));  % barrier for breaking at current y
    ebh_y(i) = ebh_func(y(i));  % barrier for healing at current y

    force(i) = inv_L(y(i)/(n+1));
    dydf(i) = ((force(i))^(-2) - (sinh(force(i)))^(-2))*(n+1);
end
dydf(1) = (n+1)/3;      % chain compliance at y = 0

ks = (n+1)*exp(-ebs_y);
kh = exp(-ebh_y);
Pe = kh./(ks+kh);

%figure;
%plot(y/n, Pe);
%figure;
%plot(force, Pe);

rho_y = zeros(1, Ny+1);
rho_f = zeros(1, Ny+1);

for i = 2:Ny
    rho_y(i) = -(Pe(i+1)-Pe(i-1))/(y(i+1)-y(i-1));
    rho_f(i) = -(Pe(i+1)-Pe(i-1))/(force(i+1)-force(i-1));
end

%rho_f2 = rho_y.*dydf;

figure
plot(y, rho_y);
figure;
plot(force, rho_f);
%hold on;
%plot(force, rho_f2);

% calculate the mean values
y_mean = 0;
f_mean = 0;
for i = 1:Ny
        y_mean = y_mean + (y(i)*rho_y(i)+y(i+1)*rho_y(i+1))/2*dy;
        f_mean = f_mean + (force(i)*rho_f(i)+force(i+1)*rho_f(i+1))/2*(force(i+1)-force(i));
end


