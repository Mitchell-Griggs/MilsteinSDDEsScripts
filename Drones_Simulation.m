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
Z0 = [0;0];
Z = zeros(2,1,N);
X3tracker = ones(2,N);
Ztracker = ones(2,N);
Zd3tracker = ones(2,N);
for j = 1:N
    if j == 1
        Z(:,1,j) = Z0+beta*X0(:,1,:)*dt;
        Zd3tracker(:,j) = Z0;
    elseif j <= 20
        Z(:,1,j) = Z(:,1,j-1) + beta*X0(:,3,:)*dt;
        Zd3tracker(:,j) = Z0;
    else
        if t(j) <= 8.5 + 1.5
            Z(:,1,j) = Z(:,1,j-1) + beta*(X(:,3,20*floor((j-1)/20))-Z(:,1,j-1))*dt;
        else
            Z(:,1,j) = Z(:,1,j-1);
        end
        Zd3tracker(:,j) = Z(:,1,j-20);
    end
    for p = 1:nd
        if j == 1
            Zd = Z0;
            X(:,p,j) = X0(:,p,1) + 0.5*v(0,tau(p),X0(:,p,1),Z0,Zd)/norm(v(0,tau(p),X0(:,p,1),Z0,Zd))*dt + sigmaj(0,tau(p),X0(:,p,1),Z0,Zd)*0;
        else    
            if j <= 21
                Zd = Z0;
            else
                Zd = Z(:,1,j-21);
            end
            X(:,p,j) = X(:,p,j-1) + 0.5*v(t(j-1),tau(p),X(:,p,j-1),Z(:,1,j-1),Zd)/norm(v(t(j-1),tau(p),X(:,p,j-1),Z(:,1,j-1),Zd))*dt + sigmaj(t(j-1),tau(p),X(:,p,j-1),Z(:,1,j-1),Zd)*(W(:,p,j)-W(:,p,j-1));
        end
    end
    X3tracker(:,j) = X(:,3,j);
    Ztracker(:,j) = Z(:,1,j);
    % The plot:
    %pause(1)  % Use this to test and see each frame slowly.
    figure(1)
    clf;hold on
    LegendX1=plot(nan,nan,'*','Color',[0,1,1]*2/3,'linewidth',10);
    LegendX2=plot(nan,nan,'*','Color',[1,0,1]*2/3,'linewidth',10);
    LegendX3=plot(nan,nan,'*','Color',[0,0,1],'linewidth',10);LX3=plot(nan,nan,'b');
    L2=plot(nan,nan,'*','Color',[1,0,0],'linewidth',10);LZ=plot(nan,nan,'r');
    L3=plot(nan,nan,'*','Color',[0,2/3,0],'linewidth',10);L4=plot(nan,nan,'--','Color',[0,2/3,0],'linewidth',1);
    L5=plot(nan,nan,'k:','LineWidth',1);
    legend([LegendX1,LegendX2,LegendX3,LX3,L2,LZ,L3,L4,L5],{'$X_1(t)$','$X_2(t)$','$X_3(t)$','path of $X_3(t)$','$Z(t)$','path of $Z(t)$','$Z(t-1)$','field of view from $Z(t-1)$','jittering axis'},'Interpreter','latex','Location','northwest','NumColumns',1,'FontSize',16);
    
    if j <= 20
        Zd = Z0;
    else
        Zd = Z(:,1,j-20);
    end
    direction = X(:,3,j) - Zd;
    
    range = norm(direction);
    theta = atan2(direction(2),direction(1))-pi/4:0.01:atan2(direction(2),direction(1))+pi/4;
    y1circle = Zd(1) + range*cos(theta); y2circle = Zd(2) + range*sin(theta);
    plot(y1circle,y2circle,'--','Color',[0,2/3,0],'linewidth',1,'HandleVisibility','off')
    perp = [-direction(2); direction(1)] / norm(direction);
    Point1 = X(:,3,j) - 2*perp;
    Point2 = X(:,3,j) + 2*perp;
    plot([Point1(1) Point2(1)], [Point1(2) Point2(2)], 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    plot(X3tracker(1,1:j),X3tracker(2,1:j),'b','LineWidth',1,'HandleVisibility','off')
    plot(Ztracker(1,1:j),Ztracker(2,1:j),'r','LineWidth',1,'HandleVisibility','off')

    % for p = 1:nd
    %     plot(X(1,p,j),X(2,p,j),'*','Color',[0,0,1],'linewidth',10,'HandleVisibility','off')
    % end
    plot(X(1,1,j),X(2,1,j),'*','Color',[0,1,1]*2/3,'linewidth',10,'HandleVisibility','off')
    plot(X(1,2,j),X(2,2,j),'*','Color',[1,0,1]*2/3,'linewidth',10,'HandleVisibility','off')
    plot(X(1,3,j),X(2,3,j),'*','Color',[0,0,1],'linewidth',10,'HandleVisibility','off')
    plot(Z(1,1,j),Z(2,1,j),'*','Color',[1,0,0],'linewidth',10,'HandleVisibility','off')
    if j <= 20
        plot(Z(1,1,1),Z(2,1,1),'*','Color',[0,2/3,0],'linewidth',10,'HandleVisibility','off')
    else
        plot(Z(1,1,j-20),Z(2,1,j-20),'*','Color',[0,2/3,0],'linewidth',10,'HandleVisibility','off')
    end

    set(gca, 'XTick', [], 'YTick', [])
    set(gcf,'position',[200,200,800,800])
    axis([-4,7,-2,9])
    xlabel('')
    ylabel('')
    hold off
    tstr = num2str(t(j), '%5.2f');
    title(['Time $t=', tstr, '$'], 'Interpreter','latex', 'FontSize',20);

    if savefigures == 1
        if t(j) == 0.2
            title('')
            saveas(gcf,'ZExample5A.eps','epsc')
        elseif t(j) == 0.5
            legend off
            title('')
            saveas(gcf,'ZExample5B.eps','epsc')
        elseif t(j) == 1
            %legend off
            title('')
            saveas(gcf,'ZExample5C.eps','epsc')
        elseif t(j) == 8
            legend off
            title('')
            saveas(gcf,'ZExample5D.eps','epsc')
        elseif t(j) == 9
            legend off
            title('')
            saveas(gcf,'ZExample5E.eps','epsc')
        elseif t(j) == 10
            legend off
            title('')
            saveas(gcf,'ZExample5F.eps','epsc')
        end
    end
end


function Vvalue = v(t,dj,X,Z,Zd)
    if t < dj
        Vvalue = X-Z;
    else
        Vvalue = X-Zd;
    end
end

function sigmajvalue = sigmaj(t,dj,X,Z,Zd)
    V = v(t,dj,X,Z,Zd);
    sigmajvalue = norm(V)^(-1)*[-V(2);V(1)];
end
