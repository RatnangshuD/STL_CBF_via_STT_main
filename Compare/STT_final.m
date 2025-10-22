% MATLAB script for CBF on a Drone with a moving ellipsoidal safe set

clc;
clear;
% close all;
% clf;

for cas = 1:4
    tic
%% Read Trajectory Data
format long;
if cas == 1
    q1 = [-1.147222222; 0.383333333; 0; 0; 0; ... 
          -0.955555556; 0; 0; 0; 0; ... 
          -1.044444444; 0.383333333; 0; 0; 0; ... 
          -0.852777778; 0; 0; 0; 0];
    q2 = [-16.87107498; 11.41027965; -2.737624136; 0.284000431; -0.010773336; ...
          41.28110288; -29.94069644; 7.523679527; -0.787750334; 0.029544048; ...
          -16.76969384; 11.41042621; -2.737366577; 0.283931461; -0.010769025; ...
          41.47972174; -30.00984186; 7.541446098; -0.78965892; 0.029616802];
    tshift = 4;
elseif cas == 2
    q1 = [-1.147222222; 0.383333333; 0; 0; 0; ...
          -0.955555556; 0; 0; 0; 0; ...
          -1.044444444; 0.383333333; 0; 0; 0; ...
          -0.852777778; 0; 0; 0; 0];
    q2 = [-61.56294811; 36.62989802; -7.874080287; 0.733347581; -0.025013661; ...
          125.7665461; -76.88066755; 16.92523909; -1.599522784; 0.055184971; ...
          -61.364337; 36.56218968; -7.856382756; 0.731326902; -0.024929093; ...
          125.9651572; -76.94106929; 16.93914191; -1.600902371; 0.055234656];
    tshift = 5;
elseif cas == 3
    q1 = [-0.970833333; 0; 0; 0; 0; ... 
          -1.139583333; 0.3375; 0; 0; 0; ... 
          -0.860416667; 0; 0; 0; 0; ... 
          -1.029166667; 0.3375; 0; 0; 0];
    q2 = [6.609233682; -3.627037902; 0.541207119; -0.023800444; 0; ...
          -8.289913916; 3.741285107; -0.503474373; 0.021746115; 0; ...
          6.712379329; -3.623893431; 0.540821853; -0.023791001; 0; ...
          -8.18676827; 3.744056067; -0.503710234; 0.021740617; 0];
    tshift = 5;
else
    q1 = [-0.983333333; 0; 0; 0; 0; ...
          -1.133333333; 0.2375; 0; 0; 0; ...
          -0.866666667; 0; 0; 0; 0; ...
          -0.98984375; 0.246354167; 0; 0; 0];
    q2 = [7.896332465; -4.228863536; 0.625985519; -0.027079889; 0; ...
          0; -0.37707672; 0.109830688; -0.006449735; 0; ...
          7.996481275; -4.222664242; 0.625487944; -0.027096203; 0; ...
          0.10014881; -0.350493552; 0.10793998; -0.006434028; 0];
end

%% Simulation & System Parameters
dt = 0.02;      % Time step [s]
tf = 9.5;        % Final time [s]
tim = 0:dt:tf;
gamma = 15.0;   % CBF parameter (class-K function gain)

% System State (simple integrator model): x_dot = u, where u is velocity
x0 = [-1; -1]; % Initial state [x; y; z]
xSTT = zeros(2, length(tim) + 1);
xSTT(:,1) = x0;

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
    CoeffLx = qnow(1:5);
    CoeffLy = qnow(6:10);
    CoeffUx = qnow(11:15);
    CoeffUy = qnow(16:20);
    
    % --- Define Time-Varying Safe Set (Ellipsoid) ---
    t_vec = [1; t; t^2; t^3; t^4];
    GamLx = CoeffLx' * t_vec;
    GamLy = CoeffLy' * t_vec;
    GamUx = CoeffUx' * t_vec;
    GamUy = CoeffUy' * t_vec;
    
    ex = ( xSTT(1,i) - 0.5*(GamUx+GamLx) ) / ( 0.5*(GamUx-GamLx) );
    ey = ( xSTT(2,i) - 0.5*(GamUy+GamLy) ) / ( 0.5*(GamUy-GamLy) );
    
    k = 2;

    ux = -k*log( (1+ex) / (1-ex) );
    uy = -k*log( (1+ey) / (1-ey) );

    uSTT = [ux; uy];
    
    % --- Update State ---
    xSTT(:,i+1) = xSTT(:,i) + dt * uSTT;
end
fprintf('Simulation finished. Generating plots...\n');

%% Plotting at specific time instances
figure(1);
subplot(1,4,cas)
% set(gcf, 'Position', [50, 50, 1500, 1000]);
hold on;

rectangle('Position', [0 -1 0.8 0.5], 'FaceColor', [1,0,0,0.5],'EdgeColor','none', FaceAlpha=0.5) % Red
rectangle('Position', [0.2 0.8 0.4 0.4], 'FaceColor', [0.4660 0.6740 0.1880],'EdgeColor','none', FaceAlpha=0.5) % Green
rectangle('Position', [-0.4 -0.4 0.8 0.8], 'FaceColor', [0,0,1,0.5],'EdgeColor','none', FaceAlpha=0.5) % Blue
rectangle('Position', [-1 -0.2 0.3 0.7], 'FaceColor', [0.9290 0.6940 0.1250],'EdgeColor','none', FaceAlpha=0.5) % Yellow


plot(xSTT(1,:), xSTT(2,:), 'k--', 'LineWidth', 2);

xlim([-1.2 1.2]); ylim([-1.2 1.2]);
grid on;
xlabel('$x_1$ (m)','Interpreter','latex','FontSize',14);
ylabel('$x_2$ (m)','Interpreter','latex','FontSize',14);
set(gca,'FontSize',16);
axis square;

fprintf('Plotting complete.\n');
toc
end