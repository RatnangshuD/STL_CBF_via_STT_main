%% cbf + 2d

clc; clear; clf;

for cas = [2,4]
tic
    % Simulation parameters
dt = 0.05; % Time step
T_sim = 5; % Simulation time
tim = 0:dt:T_sim;

% Initial condition
x0 = [-1; -1]; % Initial state
x = x0;
x_traj = x; % Store trajectory
u_traj = [];
rho_traj = [];

% Regions
T = [1; 1];
R = [0.5; -0.75];
G = [0.5; 1];
Y = [-0.75; 0.3];
dT = 0.02;
d = 0.2;

% Control parameters
alpha = @(b) b; % Class-K function
u_max = 10^7*ones(2,1); % Control saturation

% Define TV-CBF
syms x1 x2 t
r = 0.02;
gam1_0 = d - norm(x0-R) - d;
gam1_inf = (r + d)/2;
l1 = -log((r-gam1_inf) / (gam1_0 - gam1_inf));
gam1 = (gam1_0 - gam1_inf)*exp(-l1*t) + gam1_inf;

gam2_0 = d - norm(x0-G) - d;
gam2_inf = (r + d)/2;
l2 = -log((r-gam2_inf) / (gam2_0 - gam2_inf));
gam2 = (gam2_0 - gam2_inf)*exp(-l2*t) + gam2_inf;

gam3_0 = d - norm(x0-T) - d;
gam3_inf = (r + d)/2;
l3 = -log((r-gam3_inf) / (gam3_0 - gam3_inf));
gam3 = (gam3_0 - gam3_inf)*exp(-l3*t) + gam3_inf;

gam4_0 = d - norm(x0-Y) - d;
gam4_inf = (r + d)/2;
l4 = -log((r-gam4_inf) / (gam4_0 - gam4_inf));
gam4 = (gam4_0 - gam4_inf)*exp(-l4*t) + gam4_inf;

b1_sym = -gam1 + d - sqrt((x1-R(1))^2 + (x2-R(2))^2);
b2_sym = -gam2 + d - sqrt((x1-G(1))^2 + (x2-G(2))^2);
b3_sym = -gam3 + dT - sqrt((x1-T(1))^2 + (x2-T(2))^2);
b4_sym = -gam4 + dT - sqrt((x1-Y(1))^2 + (x2-Y(2))^2);

% Simulate system
for ti = tim
    if cas == 2
        if ti>=0 && ti<=2        
            b_sym = b1_sym;
        elseif ti>2 && ti<=4
            b_sym = b2_sym;
        else
            b_sym = b3_sym;
        end
    else
        if ti>=0 && ti<=2        
            b_sym = b4_sym;
        else
            b_sym = b3_sym;
        end
    end
    % Compute barrier function value at current state
    b_val = double(subs(b_sym, {x1, x2, t}, {x(1), x(2), ti}));

    % System Dynamics
    % Compute f and g matrices
    % f = [(J2-J3)/J1*x(2)*x(3); (J3-J1)/J2*x(1)*x(3); (J1-J2)/J3*x(2)*x(1)];
    % g = [1/J1 0 0;...
    %      0 1/J2 0;...
    %      0 0 1/J3];
    f = zeros(2,1);
    g = eye(2);

    % Compute constraints for QP
    A_qp = double([-subs(diff(b_sym, x1), {x1, x2, t}, {x(1), x(2), ti}), ...
                   -subs(diff(b_sym, x2), {x1, x2, t}, {x(1), x(2), ti})])*g;

    b_qp = double(-subs(diff(b_sym, t), {x1, x2, t}, {x(1), x(2), ti}) + alpha(b_val)) ...
        + double([subs(diff(b_sym, x1), {x1, x2, t}, {x(1), x(2), ti}), ...
                  subs(diff(b_sym, x2), {x1, x2, t}, {x(1), x(2), ti})])*f;

    % Solve QP for control input
    u = quadprog(eye(2), [], A_qp, b_qp, [], [], -u_max, u_max);

    % Update state
    x = x + u * dt;
    x_traj = [x_traj, x];
    u_traj = [u_traj, u];

    rho = max(d - norm(x-T), d - norm(x-G));
    rho_traj = [rho_traj, rho];
end

%% Plot results
figure(1); 

subplot(1,4,cas)
hold on; grid on;
rectangle('Position', [0 -1 0.8 0.5], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Red
rectangle('Position', [0.2 0.8 0.4 0.4], 'FaceColor', [0.4660 0.6740 0.1880],'EdgeColor','none', FaceAlpha=0.5) % Green
rectangle('Position', [-0.4 -0.4 0.8 0.8], 'FaceColor', [0,0,1,0.5],'EdgeColor','none', FaceAlpha=0.5) % Blue
rectangle('Position', [-1 -0.2 0.3 0.7], 'FaceColor', [0.9290 0.6940 0.1250],'EdgeColor','none', FaceAlpha=0.5) % Yellow

plot(x_traj(1,:), x_traj(2,:),'k-','Linewidth',1.5);
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
% hold on; grid on;
% plot(tim, x_traj(1,1:end-1), 'r', 'LineWidth', 2);
% plot(tim, x_traj(2,1:end-1), 'g', 'LineWidth', 2);
% title('State Trajectory');
% axis square 

% subplot(1,3,3)
% hold on; grid on;
% plot(tim, u_traj(1,:), 'r', 'LineWidth', 2);
% plot(tim, u_traj(2,:), 'g', 'LineWidth', 2);
% plot(tim, u_traj(3,:), 'b', 'LineWidth', 2);
% title('Control Input');
% axis square 

% subplot(1,3,3)
% hold on; grid on;
% plot(tim, rho_traj, 'r', 'LineWidth', 2);
% title('Robustness');
% axis square 

toc
end