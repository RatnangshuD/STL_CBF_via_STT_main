% CBF via STT on a Omnidirectional Mobile Robot

clc;
clear;
clf;

Coeffx = [1;0.125;0.5;-0.5;0.15625;-0.0139];
Coeffy = [1;1;-1;0.5;-0.0625;0.00097];
a1 = 0.5;
a2 = 0.25;

% State Space 
sx = linspace(-1,7);
sy = linspace(-1,7);
[sX,sY] = meshgrid(sx,sy);

% Task
dt = 1e-2;
tf = 5.92;
tim=0:dt:tf;
tplot = [0; 2; 4; 5.91];
iplot = 2;

% System
x0 = [0;0]; % initial state
xCBF = zeros(2,length(tim));
xCBF(:,1) = x0;

% Video    
% video2D = VideoWriter('test.mp4', 'MPEG-4');
% video2D.FrameRate = 60;
% video2D.Quality = 10;
% open(video2D);

%%
figure(1)
hold on;
for i=1:length(tim)
    Aplot = [0.4 0; 0 0.4];
    A = [0.3 0; 0 0.3];
    t = tim(i);
    cenx = [0, t, t^2, t^3, t^4, t^5]*Coeffx;
    ceny = [0, t, t^2, t^3, t^4, t^5]*Coeffy;
    cen = [cenx; ceny];
    
    % CBF
    dx = sX - cenx;
    dy = sY - ceny;
    Qinv = inv(A * A');
    Qinvplot = inv(Aplot * Aplot');
    
    % Evaluate h(x) = 1 - (x - c)^T Qinv (x - c) over the meshgrid
    H = 1 - (Qinv(1,1)*dx.^2 + 2*Qinv(1,2).*dx.*dy + Qinv(2,2)*dy.^2);
    Hplot = 1 - (Qinvplot(1,1)*dx.^2 + 2*Qinvplot(1,2).*dx.*dy + Qinvplot(2,2)*dy.^2);

    % CBF Controller
    gamma = 20.0;
    dX = xCBF(:,i) - cen;
    h = 1 - dX' * Qinv * dX;
    grad_h = -2 * Qinv * dX;
    
    % Set up QP
    H = eye(2);  % minimize u'u
    f = zeros(2,1);
    A_qp = -grad_h';  % for ≥ constraint
    b_qp = gamma * h / dt;
    
    % Solve QP
    options = optimoptions('quadprog', 'Display', 'off');
    uCBF = quadprog(H, f, A_qp, b_qp, [], [], [], [], [], options);

    % Plotting
    xCBF(:,i+1) = xCBF(:,i) + dt*uCBF;

    if any(t == tplot)
        figure(1)
        subplot(1,length(tplot)+1,iplot)
        contourf(sX, sY, Hplot, 100, 'LineStyle','none');
        hold on;
        rectangle("Position",[2.25,2.25,1.5,1.5],"FaceColor",'r',"FaceAlpha",0.8);
        rectangle("Position",[5,4.5,1,1],"FaceColor",'g',"FaceAlpha",0.8);
        % colorbar
        clim([-250 0]);

        plot(xCBF(1,1:i),xCBF(2,1:i),'k-','LineWidth',2); 
        plot(xCBF(1,i),xCBF(2,i),'ko'); 
        xlim([-1 7])
        ylim([-1 7])
        grid on;
        xlabel('$x_1$ (m)','Interpreter','latex','FontSize',18)
        ylabel('$x_2$ (m)','Interpreter','latex','FontSize',18)
        set(gca,'FontSize',16)
        % axis equal
        axis square
        if iplot == length(tplot)+1
            title(['$t = $', num2str(6) , ' s'],'Interpreter','latex','FontSize',18)
        else
            title(['$t = $', num2str(t) , ' s'],'Interpreter','latex','FontSize',18)
        end
        iplot = iplot+1;
    end
    
    %% 3D plot
    if mod(i,5) == 0
        figure(1)
        subplot(1,length(tplot)+1,1)
        % Generate circle points
        thetaSTT = linspace(0, 2*pi, 200);
        xSTT = cenx + 0.4*cos(thetaSTT);
        ySTT = ceny + 0.4*sin(thetaSTT);
        tSTT = t*ones(size(thetaSTT));
        fill3(xSTT, ySTT, tSTT, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 1.5);
        hold on;
        grid on;
        xlabel('$x_1$ (m)','Interpreter','latex','FontSize',18)
        ylabel('$x_2$ (m)','Interpreter','latex','FontSize',18)
        zlabel('$t$ (s)','Interpreter','latex','FontSize',18)
        set(gca,'FontSize',16)
        view([65,15]);
    end
    
    % --- End of new animation block ---

    % frame2D = getframe(1);
    % writeVideo(video2D, frame2D);
    % pause(0.01)
end
% close(video2D);