%% milp + 2d

clc;
clear;
clf;
warning('off', 'all');

for cas = 1:4
tic

% Regions
dT = 0.02;
T = [1-dT 1+dT; 1-dT 1+dT];
R = [0 0.8; -1 -0.5];
G = [0.2 0.6; 0.8 1.2];
B = [-0.4 0.4; -0.4 0.4];
Y = [-1 -0.5; 0 0.8];
d = 0.1;

% Define the time horizon and discretization
T_total = 5; % Total time
dt = 0.05; % Time step
N = T_total / dt; % Number of time steps
t_len = N;

% Define the initial state
x0 = [-1; -1];

% Define the state and control variables
x = sdpvar(2, N+1); % State variables (position)
u = sdpvar(2, N); % Control variables (velocity)

% Define the constraints
constraints = [];

% Initial condition
constraints = [constraints, x(:, 1) == x0];

% System Dynamics
for k = 1:N
    constraints = [constraints, x(:, k+1) == x(:, k) + u(:, k) * dt];
end

%% STL constraints

if cas == 1 || cas == 2
    % 1. Reach region R between time 1 and 2
    T_start_R = 1 / dt;
    T_end_R = 2 / dt;
    for k = T_start_R:T_end_R
        constraints = [constraints, (x(:, k) >= R(:, 1)+d & x(:, k) <= R(:, 2)-d)];
    end
    
    % 2. Reach region G between time 3 and 4
    T_start_G = 3 / dt;
    T_end_G = 4 / dt;
    for k = T_start_G:T_end_G
        constraints = [constraints, x(:, k) >= G(:, 1)+d & x(:, k) <= G(:, 2)-d];
    end
end

if cas == 3 || cas == 4
    % 1. Reach region Y between time 1 and 2
    T_start_Y = 1 / dt;
    T_end_Y = 2 / dt;
    for k = T_start_Y:T_end_Y
        constraints = [constraints, (x(:, k) >= Y(:, 1)+d & x(:, k) <= Y(:, 2)-d)];
    end
end

% 3. Reach region T between time 4.5 and 5
T_start_T = 4.5 / dt;
T_end_T = 5 / dt;
for k = T_start_T:T_end_T
    constraints = [constraints, x(:, k) >= T(:, 1) & x(:, k) <= T(:, 2)];
end

if cas == 1
    % 4. Avoid obstacle B for the 5 time units
    for k = 1:N
        constraints = [constraints, (x(1, k) <= B(1, 1)-d) | (x(1, k) >= B(1, 2)+d) |...
                                    (x(2, k) <= B(2, 1)-d) | (x(2, k) >= B(2, 2)+d)];
    end
end

if cas == 3
    % 4. Avoid obstacle B for the 5 time units
    for k = 1:N
        constraints = [constraints, (x(1, k) <= B(1, 1)-d) | (x(1, k) >= B(1, 2)+d) |...
                                    (x(2, k) <= B(2, 1)-d) | (x(2, k) >= B(2, 2)+d)];
    end

    % 4. Avoid obstacle G for the 5 time units
    for k = 1:N
        constraints = [constraints, (x(1, k) <= G(1, 1)-d) | (x(1, k) >= G(1, 2)+d) |...
                                    (x(2, k) <= G(2, 1)-d) | (x(2, k) >= G(2, 2)+d)];
    end
end

%% Objective: Minimize control effort
% objective = sum(sum(u.^2));
objective = sum(sum(abs(u))); % Linear approximation

% Solve the MILP problem
options = sdpsettings('solver', 'intlinprog');
% options = sdpsettings('solver', 'bmibnb');
optimize(constraints, objective, options);

% Extract and display the results
x_opt = value(x);
u_opt = value(u);

% Plot the trajectory
figure(1);
subplot(1,4,cas)
hold on; grid on;
rectangle('Position', [0 -1 0.8 0.5], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Red
rectangle('Position', [0.2 0.8 0.4 0.4], 'FaceColor', [0.4660 0.6740 0.1880],'EdgeColor','none', FaceAlpha=0.5) % Green
rectangle('Position', [-0.4 -0.4 0.8 0.8], 'FaceColor', [0,0,1,0.5],'EdgeColor','none', FaceAlpha=0.5) % Blue
rectangle('Position', [-1 -0.2 0.3 0.7], 'FaceColor', [0.9290 0.6940 0.1250],'EdgeColor','none', FaceAlpha=0.5) % Yellow

plot(x_opt(1,:), x_opt(2,:),'k-','Linewidth',1.5);
xlabel('$x_1$ (m)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
ylabel('$x_2$ (m)','interpreter','Latex','Fontsize',16,'Fontweight','bold')
grid on;
box on;
xlim([-1.2,1.2])
ylim([-1.2,1.2])
% ax = gca;
% ax.FontSize = 16;
axis square

%% Path length
% x = x_opt(1,:);
% y = x_opt(2,:);
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

toc
end