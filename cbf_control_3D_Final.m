% MATLAB script for CBF on a Drone with a moving ellipsoidal safe set

clc;
clear;
% close all;
clf;

%% Read Trajectory Data
format long;
% Attempt to read the data file, handle error if it doesn't exist
try
    q = readmatrix("Drone_1.csv");
    q1 = q(1:9);
    q2 = q(10:24);
    q3 = q(25:39);
    q4 = q(40:54);
    q1 = [q1(1); q1(2); q1(3); 0; 0; ...
          q1(4); q1(5); q1(6); 0; 0; ...
          q1(7); q1(8); q1(9); 0; 0];
catch ME
    warning('Could not read "Drone_1.csv". Using placeholder data. Error: %s', ME.message);
    % Define dummy coefficients if file is not found
    q1 = zeros(15,1); q2 = zeros(15,1); q3 = zeros(15,1); q4 = zeros(15,1);
    % Example: A simple straight line path
    q2(1:5) = [1; 0.5; 0; 0; 0]; % x = 1 + 0.5t
    q2(6:10) = [1; 0.5; 0; 0; 0]; % y = 1 + 0.5t
    q2(11:15) = [1; 0.2; 0; 0; 0]; % z = 1 + 0.2t
end

%% Simulation & System Parameters
dt = 0.05;      % Time step [s]
tf = 48;        % Final time [s]
tim = 0:dt:tf;
gamma = 15.0;   % CBF parameter (class-K function gain)

% System State (simple integrator model): x_dot = u, where u is velocity
x0 = [1; 1; 1]; % Initial state [x; y; z]
xCBF = zeros(3, length(tim) + 1);
xCBF(:,1) = x0;

%% Nominal Controller
Kp = 1.5; % Proportional gain to track the center of the safe set

%% QP Setup (for quadprog)
% min 0.5*u'*H*u + f'*u subject to A*u <= b
% We are minimizing ||u - u_nom||^2
H_qp = eye(3);
f_qp = zeros(3,1);
A_qp = zeros(1, 3); % Single CBF constraint
b_qp = zeros(1, 1);
options = optimoptions('quadprog', 'Display', 'off');

%% Main Simulation Loop
fprintf('Starting simulation...\n');
for i = 1:length(tim)
    t = tim(i);
    
    % Select appropriate trajectory segment based on time
    if t < 12
        qnow = q1;
    elseif t >= 12 && t < 24
        qnow = q2;
    elseif t >= 24 && t < 36
        qnow = q3;
    else % t >= 36
        qnow = q4;
    end
    Coeffx = qnow(1:5);
    Coeffy = qnow(6:10);
    Coeffz = qnow(11:15);
    
    % --- Define Time-Varying Safe Set (Ellipsoid) ---
    t_vec = [1; t; t^2; t^3; t^4];
    cen = [Coeffx' * t_vec; Coeffy' * t_vec; Coeffz' * t_vec];
    
    t_vec_dot = [0; 1; 2*t; 3*t^2; 4*t^3];
    cen_dot = [Coeffx' * t_vec_dot; Coeffy' * t_vec_dot; Coeffz' * t_vec_dot];
    
    A = [0.5 0 0; 0 0.5 0; 0 0 0.5]; % Spherical safe zone
    Qinv = inv(A * A');
    
    % --- Nominal Controller ---
    % Steer towards the center of the safe zone
    pos_error_nominal = cen - xCBF(:,i);
    u_nom = Kp * pos_error_nominal;
    
    % --- CBF Controller ---
    % h(x,t) = 1 - (x - c(t))' * Qinv * (x - c(t))
    dX = xCBF(:,i) - cen;
    h = 1 - dX' * Qinv * dX;
    
    % Lie derivatives for h_dot = Lf_h + Lg_h * u
    grad_h = -2 * Qinv * dX;
    
    Lf_h = grad_h' * (-cen_dot); % Part of derivative independent of control u
    Lg_h = grad_h';              % Part of derivative dependent on control u (since x_dot = u)
    
    % CBF constraint: Lf_h + Lg_h * u >= -gamma * h
    % For quadprog (A_qp * u <= b_qp):  -Lg_h * u <= Lf_h + gamma * h
    A_qp(1,:) = -Lg_h;
    b_qp(1) = Lf_h + gamma * h;
    
    % Solve QP to find safe control input uCBF closest to u_nom
    f_qp = -u_nom; % Objective: min ||u - u_nom||^2
    [uCBF, ~, exitflag] = quadprog(H_qp, f_qp, A_qp, b_qp, [], [], [], [], [], options);
    
    if exitflag ~= 1
        fprintf('QP failed at t=%.2f. Applying zero control.\n', t);
        uCBF = zeros(3,1);
    end
    
    % --- Update State ---
    xCBF(:,i+1) = xCBF(:,i) + dt * uCBF;
end
fprintf('Simulation finished. Generating plots...\n');

%% Plotting at specific time instances
figure(2);
% set(gcf, 'Position', [50, 50, 1500, 1000]);
% tplot = [0, 12, 24, 36, 48];
tplot = [4, 16, 32, 48];
dZ = 0; % Z-offset for plotting

for k = 1:length(tplot)
    subplot(1, length(tplot), k);
    hold on;
    
    % Plot static obstacles/targets
    plot_cuboid([12, 12, 0], [6, 6, 20+dZ], 'r', 0.3);
    plot_cuboid([6, 6, 6+dZ], [3, 3, 3], 'c', 0.5);
    plot_cuboid([12, 21, 9+dZ], [3, 3, 3], 'c', 0.5);
    plot_cuboid([18, 6, 8+dZ], [3, 3, 3], 'c', 0.5);
    plot_cuboid([0, 0, 0+dZ], [3, 3, 3], 'g', 0.5); % Start/Goal
    
    % Find the index corresponding to the plot time
    current_t = tplot(k);
    plot_idx = find(tim >= current_t, 1);
    if isempty(plot_idx)
        plot_idx = length(tim);
    end
    
    % Plot trajectory up to the current time
    h_traj = plot3(xCBF(1, 1:plot_idx), xCBF(2, 1:plot_idx), xCBF(3, 1:plot_idx), 'k-', 'LineWidth', 2);
    
    % Plot current drone position
    plot3(xCBF(1, plot_idx), xCBF(2, plot_idx), xCBF(3, plot_idx), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
    
    % Subplot settings
    title(['$t = $', num2str(current_t, '%.1f'), ' s'], 'Interpreter', 'latex', 'FontSize', 16);
    xlim([-5 25]); ylim([-5 25]); zlim([0 25]);
    grid on;
    xlabel('$x_1$ (m)','Interpreter','latex','FontSize',14);
    ylabel('$x_2$ (m)','Interpreter','latex','FontSize',14);
    zlabel('$x_3$ (m)','Interpreter','latex','FontSize',14);
    set(gca,'FontSize',16);
    axis square;
    view([-10,5]);
    
    % Add legend to the first subplot only
    % if k == 1
    %     legend(h_traj, 'Drone Trajectory', 'Location', 'best');
    % end
    
    hold off;
end
% sgtitle('Drone Trajectory at Different Time Instances', 'FontSize', 18, 'FontWeight', 'bold');

fprintf('Plotting complete.\n');

function plot_cuboid(origin, dims, color, alphaVal)
    % Plots a translucent cuboid in 3D space.
    x = origin(1); y = origin(2); z = origin(3);
    dx = dims(1); dy = dims(2); dz = dims(3);
    
    vertices = [x,y,z; x+dx,y,z; x+dx,y+dy,z; x,y+dy,z; ...
                x,y,z+dz; x+dx,y,z+dz; x+dx,y+dy,z+dz; x,y+dy,z+dz];
    
    faces = [1,2,3,4; 5,6,7,8; 1,2,6,5; 2,3,7,6; 3,4,8,7; 4,1,5,8];
    
    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', color, 'FaceAlpha', alphaVal, ...
          'EdgeColor', 'k', 'LineWidth', 0.5);
end

