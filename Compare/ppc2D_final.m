%% ppc + 2d

clc;
clear;
clf;
global cas
for cas = [2,4]

X_init = [-1; -1];
t_len = 500;
t_span = linspace(0,5,t_len);
tic
[t, X] = ode45(@omni, t_span, X_init);
toc

%% plot
figure(1)
hold on;
subplot(1,4,cas)
hold on; grid on;
rectangle('Position', [0 -1 0.8 0.5], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Red
rectangle('Position', [0.2 0.8 0.4 0.4], 'FaceColor', [0.4660 0.6740 0.1880],'EdgeColor','none', FaceAlpha=0.5) % Green
rectangle('Position', [-0.4 -0.4 0.8 0.8], 'FaceColor', [0,0,1,0.5],'EdgeColor','none', FaceAlpha=0.5) % Blue
rectangle('Position', [-1 -0.2 0.3 0.7], 'FaceColor', [0.9290 0.6940 0.1250],'EdgeColor','none', FaceAlpha=0.5) % Yellow

plot(X(:,1), X(:,2),'k-','Linewidth',1.5);
xlabel('$x_1$ (m)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
ylabel('$x_2$ (m)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
grid on;
box on;
xlim([-1.2,1.2])
ylim([-1.2,1.2])
% ax = gca;
% ax.FontSize = 16;
axis square

% subplot(1,3,2)
% plot(t,X, 'Linewidth', 1.5)
% grid on;
% legend('x', 'y')

end
%% Path length
% x = X(:,1);
% y = X(:,2);
% 
% % Smooth the data using Savitzky-Golay filter
% window_size = 5; %round(t_len/200)*2+1; % Choose an appropriate window size
% order = 3;        % Polynomial order
% x_smooth = sgolayfilt(x, order, window_size);
% y_smooth = sgolayfilt(y, order, window_size);
% 
% % Compute path length
% dx = diff(x_smooth);
% dy = diff(y_smooth);
% segment_lengths = sqrt(dx.^2 + dy.^2);
% path_length = sum(segment_lengths);
% 
% % Display results
% fprintf('Smoothed Path Length: %.4f\n', path_length);
% 
% toc

%% System Dynamics
function dXdt = omni(t,X)
u = real(control(t,X));
dXdt = u;
end

%% Control Law
function u = control(t,X)
global cas

T = [1; 1];
R = [0.5; -0.75];
G = [0.5; 1];
Y = [-0.75; 0.3];
dT = 0.02;
d = 0.2;

x = X(1);
y = X(2);

if cas == 2
    if t>=0 && t<=2
        rho = d - sqrt((x-R(1))^2 + (y-R(2))^2);
        drho_dx = -[x-R(1); y-R(2)]/sqrt((x-R(1))^2 + (y-R(2))^2);
        spec = 1;
    elseif t>2 && t<=4
        rho = d - sqrt((x-G(1))^2 + (y-G(2))^2);
        drho_dx = -[x-G(1); y-G(2)]/sqrt((x-G(1))^2 + (y-G(2))^2);
        spec = 2;
    else
        rho = dT - sqrt((x-T(1))^2 + (y-T(2))^2);
        drho_dx = -[x-T(1); y-T(2)]/sqrt((x-T(1))^2 + (y-T(2))^2);
        spec = 3;
    end
else
    if t>=0 && t<=2
        rho = d - sqrt((x-Y(1))^2 + (y-Y(2))^2);
        drho_dx = -[x-Y(1); y-Y(2)]/sqrt((x-Y(1))^2 + (y-Y(2))^2);
        spec = 1;
    % elseif t>2 && t<=4
    %     rho = d - sqrt((x-G(1))^2 + (y-G(2))^2);
    %     drho_dx = -[x-G(1); y-G(2)]/sqrt((x-G(1))^2 + (y-G(2))^2);
    %     spec = 2;
    else
        rho = dT - sqrt((x-T(1))^2 + (y-T(2))^2);
        drho_dx = -[x-T(1); y-T(2)]/sqrt((x-T(1))^2 + (y-T(2))^2);
        spec = 3;
    end
end

eps = err(t,rho,spec,cas);

u = -10*eps*drho_dx;
end


%% error
function eps = err(t,rho,spec,cas)

%% Funnel
gamma_0 = 1;
gamma_inf = 0.1;
rho_max = 0.5;

t_f1 = 0; % funnel 1 start time
t_f2 = 2; % funnel 2 start time
t_f3 = 4; % funnel 3 start time

t1_star = 1;
l1 = -1/t1_star*log((rho_max-gamma_inf)/(gamma_0-gamma_inf));
gamma1 = (gamma_0-gamma_inf)*exp(-l1*(t-t_f1))+gamma_inf;

t2_star = 3;
l2 = -1/t2_star*log((rho_max-gamma_inf)/(gamma_0-gamma_inf));
gamma2 = (gamma_0-gamma_inf)*exp(-l2*(t-t_f2))+gamma_inf;

t3_star = 4.5;
l3 = -1/t3_star*log((rho_max-gamma_inf)/(gamma_0-gamma_inf));
gamma3 = (gamma_0-gamma_inf)*exp(-l3*(t-t_f3))+gamma_inf;

%%
if spec == 1
    e = (rho-rho_max)/gamma1;
elseif spec == 2
    e = (rho-rho_max)/gamma2;
else
    e = (rho-rho_max)/gamma3;
end
eps = log(-(1+e)/(e));
end

%% 2-Norm
function f = normTwo(a,b) 
    x = a(1)-b(1);
    y = a(2)-b(2);
    f = sqrt(x^2+y^2);
end

