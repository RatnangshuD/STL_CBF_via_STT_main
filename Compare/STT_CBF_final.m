% MATLAB script for CBF on a Drone with a moving ellipsoidal safe set

clc;
clear;
% close all;
% clf;

for cas = 1:4
    tic
%% Read Trajectory Data
format long;
tshift = 5;
if cas == 1
    q1 = [-1.04583; 0.383333; 0; 0; 0; ...
          -0.95417; 0.019444; 0; 0; 0];
    q2 = [-2.198685043; 1.155933258; -0.123150965; -0.004639819; 0.000847972; ...
       11.0383636; -10.48765614; 3.143331109; -0.369858362; 0.015123053];
elseif cas == 2
    q1 = [-1.01388888888888; 0.255555555555555; 0; 0; 0; ...          
          -0.986111111111111; 0.0271604938271604; 0; 0; 0];
    q2 = [-18.49604552; 9.028360768; -1.411263289; 0.071197702; 0; ...
          8.975683728; -4.788092516; 0.730685823; -0.033221403; 0];
elseif cas == 3
    q1 = [-0.970833333; 0.030555556; 0; 0; 0; ...
          -1.029166667; 0.316666667; 0; 0; 0];
    q2 = [5.872748112; -3.241756562; 0.487695074; -0.021395182; 0; ...
          -10.70864347; 5.043871359; -0.72934466; 0.034216559; 0];
else
    q1 = [-0.9625; 0.0166666666666666; 0; 0; 0; ...
          -0.9625; 0.194444444444444; 0; 0; 0];
    q2 = [1.20313786; -0.849588477; 0.086625514; 0; 0; ...
          -0.881910151; 0.162208505; 0.003223594; 0; 0];
end

%% Simulation & System Parameters
dt = 0.02;      % Time step [s]
tf = 9.5;        % Final time [s]
tim = 0:dt:tf;
gamma = 15.0;   % CBF parameter (class-K function gain)

% System State (simple integrator model): x_dot = u, where u is velocity
x0 = [-1; -1]; % Initial state [x; y; z]
xCBF = zeros(2, length(tim) + 1);
xCBF(:,1) = x0;

%% Nominal Controller
Kp = 1.5; % Proportional gain to track the center of the safe set

%% QP Setup (for quadprog)
% min 0.5*u'*H*u + f'*u subject to A*u <= b
% We are minimizing ||u - u_nom||^2
H_qp = eye(2);
f_qp = zeros(2,1);
A_qp = zeros(1, 2); % Single CBF constraint
b_qp = zeros(1, 1);
options = optimoptions('quadprog', 'Display', 'off');

%% Main Simulation Loop
fprintf('Starting simulation...\n');
for i = 1:length(tim)
    t = tim(i);
    
    % Select appropriate trajectory segment based on time
    if t < tshift
        qnow = q1;
    elseif t >= tshift
        qnow = q2;
    else % t >= tshift
        qnow = q2;
    end
    Coeffx = qnow(1:5);
    Coeffy = qnow(6:10);
    
    % --- Define Time-Varying Safe Set (Ellipsoid) ---
    t_vec = [1; t; t^2; t^3; t^4];
    cen = [Coeffx' * t_vec; Coeffy' * t_vec];
    
    t_vec_dot = [0; 1; 2*t; 3*t^2; 4*t^3];
    cen_dot = [Coeffx' * t_vec_dot; Coeffy' * t_vec_dot];
    
    A = [0.1 0; 0 0.1]; % Spherical safe zone
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
figure(1);
subplot(1,4,cas)
% set(gcf, 'Position', [50, 50, 1500, 1000]);
hold on;

% rectangle('Position', [0 -1 0.8 0.5], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Red
% rectangle('Position', [0.2 0.8 0.4 0.4], 'FaceColor', [0.4660 0.6740 0.1880],'EdgeColor','none', FaceAlpha=0.5) % Green
% rectangle('Position', [-0.4 -0.4 0.8 0.8], 'FaceColor', [0,0,1,0.5],'EdgeColor','none', FaceAlpha=0.5) % Blue
% rectangle('Position', [-1 -0.2 0.3 0.7], 'FaceColor', [0.9290 0.6940 0.1250],'EdgeColor','none', FaceAlpha=0.5) % Yellow


plot(xCBF(1,:), xCBF(2,:), 'k-', 'LineWidth', 2);

xlim([-1.2 1.2]); ylim([-1.2 1.2]);
grid on;
xlabel('$x_1$ (m)','Interpreter','latex','FontSize',14);
ylabel('$x_2$ (m)','Interpreter','latex','FontSize',14);
set(gca,'FontSize',16);
axis square;

fprintf('Plotting complete.\n');
toc
end