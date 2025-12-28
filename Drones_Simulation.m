savefigures = 0;  % Set this to 0 to avoid saving figures, which saves significant computation time.
rng('default')
range = 2;
T = 10;
dt = 0.05;
N = T/dt;
t = dt:dt:T;
nd = 3; % Number of drones.
tau = (1:nd)/nd;
W = cumsum(sqrt(dt)*randn(1,nd,N),3);
beta = 1/1;

% Drone movements:
X0 = 10*rand(2,nd,1)-5;
X0(:,1,1) = [1;1]/3;
X0(:,2,1) = [0;1]/3;
X0(:,3,1) = [1;0]/5;
X = zeros(2,nd,N);
Y0 = [0;0];
Y = zeros(2,1,N);
X3tracker = ones(2,N);
Ytracker = ones(2,N);
Yd3tracker = ones(2,N);
for j = 1:N
    if j == 1
        Y(:,1,j) = Y0+beta*X0(:,1,:)*dt;
        Yd3tracker(:,j) = Y0;
    elseif j <= 20
        Y(:,1,j) = Y(:,1,j-1) + beta*X0(:,3,:)*dt;
        Yd3tracker(:,j) = Y0;
    else
        if t(j) <= 8.5
            Y(:,1,j) = Y(:,1,j-1) + beta*(X(:,3,20*floor((j-1)/20))-Y(:,1,j-1))*dt;
        else
            Y(:,1,j) = Y(:,1,j-1);
        end
        Yd3tracker(:,j) = Y(:,1,j-20);
    end
    for p = 1:nd
        if j == 1
            Yd = Y0;
            X(:,p,j) = X0(:,p,1) + 0.5*v(0,tau(p),X0(:,p,1),Y0,Yd)/norm(v(0,tau(p),X0(:,p,1),Y0,Yd))*dt + sigmaj(0,tau(p),X0(:,p,1),Y0,Yd)*0;
        else    
            if j <= 21
                Yd = Y0;
            else
                Yd = Y(:,1,j-21);
            end
            X(:,p,j) = X(:,p,j-1) + 0.5*v(t(j-1),tau(p),X(:,p,j-1),Y(:,1,j-1),Yd)/norm(v(t(j-1),tau(p),X(:,p,j-1),Y(:,1,j-1),Yd))*dt + sigmaj(t(j-1),tau(p),X(:,p,j-1),Y(:,1,j-1),Yd)*(W(:,p,j)-W(:,p,j-1));
        end
    end
    X3tracker(:,j) = X(:,3,j);
    Ytracker(:,j) = Y(:,1,j);
    % The plot:
    figure(1)
    clf;hold on
    L1=plot(nan,nan,'*','Color',[0,0,1],'linewidth',10);LX3=plot(nan,nan,'b');
    L2=plot(nan,nan,'*','Color',[1,0,0],'linewidth',10);LY=plot(nan,nan,'r');
    L3=plot(nan,nan,'*','Color',[0,2/3,0],'linewidth',10);L4=plot(nan,nan,'--','Color',[0,2/3,0],'linewidth',1);
    L5=plot(nan,nan,'k:','LineWidth',1);
    legend([L1,LX3,L2,LY,L3,L4,L5],{'$X_j(t)$, $j=1,2,3$','$Y(t)$','$Y(t-1)$','$X_3(t)$ Path','$Y(t)$ Path','$Y(t-1)$ Field of View','Jittering Axis'},'Interpreter','latex','Location','northwest','NumColumns',1,'FontSize',16);
    
    if j <= 20
        Yd = Y0;
    else
        Yd = Y(:,1,j-20);
    end
    direction = X(:,3,j) - Yd;
    
    range = norm(direction);
    theta = atan2(direction(2),direction(1))-pi/4:0.01:atan2(direction(2),direction(1))+pi/4;
    y1circle = Yd(1) + range*cos(theta); y2circle = Yd(2) + range*sin(theta);
    plot(y1circle,y2circle,'--','Color',[0,2/3,0],'linewidth',1,'HandleVisibility','off')
    perp = [-direction(2); direction(1)] / norm(direction);
    Point1 = X(:,3,j) - 2*perp;
    Point2 = X(:,3,j) + 2*perp;
    plot([Point1(1) Point2(1)], [Point1(2) Point2(2)], 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    plot(X3tracker(1,1:j),X3tracker(2,1:j),'b','LineWidth',1,'HandleVisibility','off')
    plot(Ytracker(1,1:j),Ytracker(2,1:j),'r','LineWidth',1,'HandleVisibility','off')

    for p = 1:nd
        plot(X(1,p,j),X(2,p,j),'*','Color',[0,0,1],'linewidth',10,'HandleVisibility','off')
    end
    plot(Y(1,1,j),Y(2,1,j),'*','Color',[1,0,0],'linewidth',10,'HandleVisibility','off')
    if j <= 20
        plot(Y(1,1,1),Y(2,1,1),'*','Color',[0,2/3,0],'linewidth',10,'HandleVisibility','off')
    else
        plot(Y(1,1,j-20),Y(2,1,j-20),'*','Color',[0,2/3,0],'linewidth',10,'HandleVisibility','off')
    end

    set(gca, 'XTick', [], 'YTick', [])
    set(gcf,'position',[200,200,800,800])
    axis([-8,9,-4,13])
    xlabel('')
    ylabel('')
    hold off
    tstr = num2str(t(j), '%5.2f');
    title(['Time $t=', tstr, '$'], 'Interpreter','latex', 'FontSize',20);

    if savefigures == 1
        if t(j) == 0.2
            title('')
            saveas(gcf,'DronesFigure1.eps','epsc')
        elseif t(j) == 0.5
            legend off
            title('')
            saveas(gcf,'DronesFigure2.eps','epsc')
        elseif t(j) == 1
            legend off
            title('')
            saveas(gcf,'DronesFigure3.eps','epsc')
        elseif t(j) == 8
            legend off
            title('')
            saveas(gcf,'DronesFigure4.eps','epsc')
        elseif t(j) == 9
            legend off
            title('')
            saveas(gcf,'DronesFigure5.eps','epsc')
        elseif t(j) == 10
            legend off
            title('')
            saveas(gcf,'DronesFigure6.eps','epsc')
        end
    end
end


function Vvalue = v(t,dj,X,Y,Yd)
    if t < dj
        Vvalue = X-Y;
    else
        Vvalue = X-Yd;
    end
end

function sigmajvalue = sigmaj(t,dj,X,Y,Yd)
    V = v(t,dj,X,Y,Yd);
    sigmajvalue = norm(V)^(-1)*[-V(2);V(1)];
end