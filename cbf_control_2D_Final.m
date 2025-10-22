% MATLAB script for CBF on a Differential Drive Robot with a moving safe set


% close all;
clc;
clear;

for cas = 1:2

%% Simulation & Robot Parameters
dt = 0.02;      % Time step [s]
tf = 15;        % Final time [s]
tim = 0:dt:tf;
i0 = 10;
gamma = 10.0;   % CBF parameter

% System State: [x; y; theta]
q0 = [0; 0; pi/4]; % Initial state
qCBF = zeros(3, length(tim) + 1);
qCBF(:,1) = q0;

%% Nominal Controller Gains (to follow the center of the safe set)
Kp = 1.0; % Proportional gain for linear velocity
Ka = 2.0; % Proportional gain for angular velocity

%% Trajectory Coefficients for the Center of the Safe Set
% These define the path c(t) = [cenx(t); ceny(t)]
if cas == 1
    clf;
    % DDR_T1
    Coeffx = [1; -0.00293; 0.200195; -0.00781; 0];
    Coeffy = [1; -0.9685; 0.5; -0.03363; 0.000629];
else
    % DDR_T2
    Coeffx = [1; 2.625; -0.09375; 0; 0];
    Coeffy = [1; 0.75; 0; 0.007813; -0.00043];
end
%% Plotting Setup
figure(1)
% set(gcf, 'Position', [100, 100, 1600, 600]); % Set figure size
% tplot = [0, 3, 6, 9, 12, 15]; % Timestamps for 2D plots
tplot = 15; % Timestamps for 2D plots
iplot = 2; % Subplot index counter

% State Space for contour plots
sx = linspace(-1, 22);
sy = linspace(-1, 22);
[sX, sY] = meshgrid(sx, sy);

%% QP Setup (for quadprog)
% min 0.5*u'*H*u + f'*u subject to A*u <= b
% We are minimizing ||u - u_nom||^2
H_qp = eye(2);
f_qp = zeros(2,1);

% Define constraints: only the single CBF constraint
A_qp = zeros(1, 2);
b_qp = zeros(1, 1);
options = optimoptions('quadprog', 'Display', 'off');

%% Main Simulation Loop
fprintf('Starting simulation...\n');
for i = 1:length(tim)
    t = tim(i);
    
    % Current robot state
    q = qCBF(:,i);
    x_pos = q(1);
    y_pos = q(2);
    theta = q(3);
    
    % --- Define the Time-Varying Safe Set ---
    t_vec = [1, t, t^2, t^3, t^4]'; % Use 1 for constant term
    cen = [Coeffx' * t_vec; Coeffy' * t_vec];
    
    % Time derivative of the center, c_dot
    t_vec_dot = [0, 1, 2*t, 3*t^2, 4*t^3]';
    cen_dot = [Coeffx' * t_vec_dot; Coeffy' * t_vec_dot];
    
    % Ellipsoid shape matrix for h(x) = 1 - (x-c)'*Qinv*(x-c)
    A = [0.5 0; 0 0.5];
    Qinv = inv(A * A');

    % --- Nominal Controller ---
    % A simple controller to track the center of the safe set 'cen'
    dist_error = norm([x_pos; y_pos] - cen);
    angle_to_cen = atan2(cen(2) - y_pos, cen(1) - x_pos);
    angle_error = wrapToPi(angle_to_cen - theta);

    v_nom = Kp * dist_error;
    omega_nom = Ka * angle_error;
    u_nom = [v_nom; omega_nom];
    
    % --- CBF Controller ---
    % h(q,t) = 1 - (pos - cen)' * Qinv * (pos - cen)
    pos_error = [x_pos; y_pos] - cen;
    h = 1 - pos_error' * Qinv * pos_error;
    
    % Lie derivatives for differential drive model
    % dh/dt = (dh/dx)*x_dot + dh/dt
    % dh/dt = Lf_h + Lg_h * u
    grad_h_pos = -2 * Qinv * pos_error; % Gradient of h w.r.t position [x,y]
    
    % Part of derivative not dependent on control u
    Lf_h = grad_h_pos' * (-cen_dot);
    
    % Part of derivative dependent on control u = [v; omega]
    % Lg_h = [dh/dx * [cos(theta); sin(theta)], 0]
    Lg_h = [grad_h_pos(1)*cos(theta) + grad_h_pos(2)*sin(theta), 0];

    % CBF constraint: Lf_h + Lg_h * u >= -gamma * h
    % For quadprog (A*u <= b):  -Lg_h * u <= gamma*h + Lf_h
    A_qp(1,:) = -Lg_h;
    b_qp(1) = gamma * h + Lf_h;
    
    % Solve QP
    f_qp = -u_nom; % To minimize ||u-u_nom||^2
    [uCBF, ~, exitflag] = quadprog(H_qp, f_qp, A_qp, b_qp, [], [], [], [], [], options);
    
    if exitflag ~= 1
        fprintf('QP failed at t=%.2f. Applying zero control.\n', t);
        uCBF = [0; 0];
    end
    
    % --- Update State using Differential Drive Kinematics ---
    v = uCBF(1);
    omega = uCBF(2);
    q_dot = [v * cos(theta);
             v * sin(theta);
             omega];
    
    qCBF(:, i+1) = q + dt * q_dot;
    qCBF(3, i+1) = wrapToPi(qCBF(3, i+1)); % Keep theta in [-pi, pi]
    
    % --- Plotting ---
    % Check if current time is in the list of times to plot
    if any(abs(t - tplot) < dt/2)
        % subplot(2, length(tplot)/2 + 1, iplot);
        subplot(1, 4, iplot + (cas-1)*2);
        
        % Calculate h over the meshgrid for visualization
        dx = sX - cen(1);
        dy = sY - cen(2);
        H_grid = 1 - (Qinv(1,1)*dx.^2 + 2*Qinv(1,2).*dx.*dy + Qinv(2,2)*dy.^2);
        
        % contourf(sX, sY, H_grid, 100, 'LineStyle','none');
        clim([-2000 0]);
        hold on;
        
        % Plot the zero-level set of h (the boundary)
        % contour(sX, sY, H_grid, [0 0], 'r-', 'LineWidth', 2);
        
        % Specifications
        rectangle("Position",[0,0,3,3],"FaceColor",'c',"FaceAlpha",0.4); % Start
        rectangle("Position",[9,6,3,3],"FaceColor",'r',"FaceAlpha",0.8); % Obstacle
        rectangle("Position",[6,6,3,3],"FaceColor",'c',"FaceAlpha",0.8); % T1
        rectangle("Position",[12,6,3,3],"FaceColor",'c',"FaceAlpha",0.8); % T2
        rectangle("Position",[18,15,3,3],"FaceColor",'g',"FaceAlpha",0.8); % G
        
        % Plot the robot's trajectory so far
        plot(qCBF(1,i0:i), qCBF(2,i0:i), 'k-', 'LineWidth', 2); 
        
        % Plot robot's current position and orientation as a triangle
        robot_size = 1.5; % Size of the triangle
        tip_x = x_pos + robot_size * cos(theta);
        tip_y = y_pos + robot_size * sin(theta);
        
        base_angle1 = theta + 3*pi/4; % Angle for one base corner
        base_angle2 = theta - 3*pi/4; % Angle for other base corner
        
        base_x1 = x_pos + robot_size/2 * cos(base_angle1);
        base_y1 = y_pos + robot_size/2 * sin(base_angle1);
        base_x2 = x_pos + robot_size/2 * cos(base_angle2);
        base_y2 = y_pos + robot_size/2 * sin(base_angle2);
        
        fill([tip_x, base_x1, base_x2], [tip_y, base_y1, base_y2], 'k', 'EdgeColor', 'k', 'LineWidth', 1);
        
        xlim([-1 22]);
        ylim([-1 22]);
        grid on;
        xlabel('$x_1$ (m)','Interpreter','latex','FontSize',14);
        ylabel('$x_2$ (m)','Interpreter','latex','FontSize',14);
        set(gca,'FontSize',18);
        axis square;
        title(['$t = $', num2str(t, '%.1f'), ' s'], 'Interpreter', 'latex', 'FontSize', 16);
        iplot = iplot + 1;
        % if iplot == 5, iplot = iplot + 1; end % Skip middle plot area
    end
    
    % 3D plot update
    if mod(i, 5) == 0 % Update less frequently to speed up
        % subplot(2, length(tplot)/2 + 1, [1, 5]); % Use the large subplot area
        subplot(1, 4, 1 + (cas-1)*2);
        
        hold on;
        
        % Plot the moving safe set boundary as a tube
        thetaSTT = linspace(0, 2*pi, 200);
        xSTT = cen(1) + 1*cos(thetaSTT);
        ySTT = cen(2) + 1*sin(thetaSTT);
        tSTT = t*ones(size(thetaSTT));
        % plot3(xSTT, ySTT, tSTT, 'b-', 'LineWidth', 2);
        fill3(xSTT, ySTT, tSTT, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'LineWidth', 1.5);
        
        % Plot trajectory in 3D (x, y, t)
        plot3(qCBF(1,i0:i), qCBF(2,i0:i), tim(i0:i), 'k-', 'LineWidth', 2);

        grid on;
        xlabel('$x_1$ (m)','Interpreter','latex','FontSize',14);
        ylabel('$x_2$ (m)','Interpreter','latex','FontSize',14);
        zlabel('$t$ (s)','Interpreter','latex','FontSize',14);
        set(gca,'FontSize',18);
        view([65, 15]);
        xlim([-1 22]); ylim([-1 22]); zlim([0 tf]);
        hold off;
    end
end
fprintf('Simulation finished.\n');

end