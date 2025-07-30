% This script is for plotting the mean-square error graphs of the semilinear,
% homogeneous, 2-dimensional SDDE
%  dX(t) = \sum_{j=0}^m ( A_j*X(t) + B_j(X(t),X(t-d1),X(t-d2)) dW_j(t),
% 0<=t<=T, with X(t) = X(0)=X0 for t<=0,   where m = 2 and d = 2, and the
% delays d1 and d2 are not necessarily on the scheme mesh, but the schemes have initial step size h0 that divides T.

%% USER PARAMETER:

runs = 5;   % Number of trajectories.

% Please feel free to change the "runs" (the number of trials) at will. For
% example, verify that this script works with "runs = 10."


%% OTHER PARAMETERS:

% SIMULATION OPTIONS:
useRefinedMilsteinschemes = 1; % Set this to 1 to simulate the refined Milstein schemes and the double integrals. This is significantly slower than without.
saveplots = 1;

% Simulation parameters:
dt = 2^-10;    % Smallest time increment in the approximations.
refdt = 2^-13; % Step size in reference solution.
runstosave = 1000;   % Number of trials in order to save the plots.

% Equation parameters:
X0 = [0.8;0.2];
A0 = [-0.1,0.03;-0.2,-0.04];   A1 = [0.05,0.04;0.02,0.03];  A2 = [0.05,0.03;0.04,0.01];

B0 = @(X,Xd1,Xd2) [sin(X(1));cos(X(2))]/5; 
B1 = @(X,Xd1,Xd2) [Xd2(1)-Xd1(1);Xd1(2)-Xd2(2)]/3; 
B2 = @(X,Xd1,Xd2) [exp(-X(2)^2)+exp(-Xd1(1)^2)+exp(-Xd1(2)^2);exp(-X(1)^2)+exp(-Xd2(1)^2)+exp(-Xd2(2)^2)]/10;

%Calculate Jacobians:
syms x1 x2 y1 y2 z1 z2
    B1xtemp = eval(['@(x1,x2)' char(jacobian(B1([x1;x2],[y1;y2],[z1;z2]),[x1;x2]))]);
    B1x = @(X,Y,Z) B1xtemp(X(1),X(2));
    B1xd1temp = eval(['@(y1,y2)' char(jacobian(B1([x1;x2],[y1;y2],[z1;z2]),[y1;y2]))]);
    B1xd1 = @(X,Y,Z) B1xd1temp(Y(1),Y(2));
    B1xd2temp = eval(['@(z1,z2)' char(jacobian(B1([x1;x2],[y1;y2],[z1;z2]),[z1;z2]))]);
    B1xd2 = @(X,Y,Z) B1xd2temp(Z(1),Z(2));
    B2xtemp = eval(['@(x1,x2)' char(jacobian(B2([x1;x2],[y1;y2],[z1;z2]),[x1;x2]))]);
    B2x = @(X,Y,Z) B2xtemp(X(1),X(2));
    B2xd1temp = eval(['@(y1,y2)' char(jacobian(B2([x1;x2],[y1;y2],[z1;z2]),[y1;y2]))]);
    B2xd1 = @(X,Y,Z) B2xd1temp(Y(1),Y(2));
    B2xd2temp = eval(['@(z1,z2)' char(jacobian(B2([x1;x2],[y1;y2],[z1;z2]),[z1;z2]))]);
    B2xd2 = @(X,Y,Z) B2xd2temp(Z(1),Z(2));
clear x1 x2 y1 y2 z1 z2

phi = @(t) [0.8+0.2.*t.^2.*sin(3*pi*t);0.2+0.2.*t.^2.*cos(2*pi*t)];

T = 4;  % T > 1 is required in this script.

%% CALCULATIONS BEGIN HERE.

simulationtime = tic;

for simulation = 1:3

if simulation == 1
    d1 = 1; d2 = 1/2;
elseif simulation == 2
    d1 = 1; d2 = pi/4;
elseif simulation == 3
    d1 = exp(2)/5; d2 = pi/4;
elseif simulation == 4
    d1 = 1/2; d2 = pi/4;
elseif simulation == 5
    d1 = 8/7; d2 = pi/4;
elseif simulation == 6
    d1 = 3/4; d2 = pi/7;
end

% Times for observations of the error graphs:
tobservations = [T];%[0.25,d1,3*d2,4];

NrValuesPotentially = 1:100; % The upper bound Nr on this needs to have Nr*r >= T.

% Observation parameters:
showerrorplots = 1;
hvalues = 2.^(min([floor(log(d1)/log(2)),floor(log(d2)/log(2)),0]):-1:log(dt)/log(2));        % Set of step sizes to look at.
samplepath = 1;                             % Sample trajectory to view.
samplecomponent = 1;                        % Sample component to view.
samplehvalue = 1/2;                          % Sample step size to view.

%% INITIALISING CALCULATIONS:

rng('default')
T = ceil(T);    % This is ensuring T is an integer.

% Calculate the initial mesh:
mesh0 = refdt:refdt:T;
% Set the observation times:
N1 = find(NrValuesPotentially*d1>=T,1);
N2 = find(NrValuesPotentially*d2>=T,1);
TO = [mesh0,(1:N1)*d1,(1:N2)*d2,tobservations];
TO = sort(TO);    TO = TO(0<TO&TO<=T);    TO = unique(TO);

% Set the complete mesh:
TC = TO;
for i = 0:N1
    for j = 0:N2
        TC = [TC,TC - i*d1 - j*d2];
        TC = sort(TC);    TC = TC(0<TC&TC<=T);    TC = unique(TC);
    end
end
TC = sort(TC);    TC = TC(0<TC&TC<=T);    TC = unique(TC);

t = TC; % Time.

% Find the index positions in t when t=k*r1 for some integer k:
td1index = zeros(1,N1);
for i = 1:N1
    if i*d1 <= T
        td1index(i) = find(t==i*d1,1);
    end
end
td1index = td1index(td1index~=0);

% Find the index positions in t when t=k*r2 for some integer k:
td2index = zeros(1,N2);
for i = 1:N2
    if i*d2 <= T
        td2index(i) = find(t==i*d2,1);
    end
end
td2index = td2index(td2index~=0);

% Find the index positions in t to observe the error graphs at:
tindexobservations = zeros(1,length(tobservations));
for i = 1:length(tobservations)
    tindexobservations(i) = find(t==tobservations(i));
end
NO = length(tindexobservations); % Number of observation times.

% To test these, run the script and then put t(tindexobservations) into the command window.


A0hat = A0 - (A1^2+A2^2)/2;

Xrefvalues = zeros(2,NO,runs);

LieA0A1 = A0*A1-A1*A0;
LieA0A2 = A0*A2-A2*A0;
LieA1A2 = A1*A2-A2*A1;

% Initialise the errors, and compile them in vectors, so that we can view
% their means plotted against the step sizes, at the end of this script.
% The first variable refers to the time that the error is observed at. For
% example, ErrorEM(1,.,.) is the error at t = 1, while ErrorEM(2,.,.) is the
% error at t = 2.
ErrorEM = ones(NO,runs,length(hvalues));     % Strong errors compared with X, a reference solution formed by the Milstein scheme with varying step size, defined in Chapter 6 of my thesis.
ErrorMilSim = ones(NO,runs,length(hvalues));
ErrorMEM = ones(NO,runs,length(hvalues));
ErrorMMSim = ones(NO,runs,length(hvalues));
ErrorMilRef = ones(NO,runs,length(hvalues));    
ErrorMMRef = ones(NO,runs,length(hvalues));

% Stacks for the values:
EMvalues = zeros(2,NO,runs,length(hvalues));  % Indexing is time, then trial, then step size for the scheme.
MilSimvalues = zeros(2,NO,runs,length(hvalues));
MEMvalues = zeros(2,NO,runs,length(hvalues));
MMSimvalues = zeros(2,NO,runs,length(hvalues));
MilRefvalues = zeros(2,NO,runs,length(hvalues));
MMRefvalues = zeros(2,NO,runs,length(hvalues));

MSErrorEM = ones(NO,length(hvalues));
MSErrorMilSim = ones(NO,length(hvalues));
MSErrorMEM = ones(NO,length(hvalues));
MSErrorMMSim = ones(NO,length(hvalues));
MSErrorMilRef = ones(NO,length(hvalues));
MSErrorMMRef = ones(NO,length(hvalues));

%% SIMULATIONS BEGIN HERE:

%% SIMULATIONS BEGIN HERE:

for trial = 1:runs

    if mod(trial,100)==0; disp(['Trial number: ',num2str(trial),'     Simulation number: ',num2str(simulation)]); end

%% SIMULATION OF WIENER PATH:

    % W1 and W1d, beginning with W1:

    % Initialise:
    dW1initial = sqrt(t(1))*randn(1,1);
    dW2initial = sqrt(t(1))*randn(1,1);
    dW1 = zeros(1,length(t));
    dW2 = zeros(1,length(t));
    W1 = zeros(1,length(t));
    W2 = zeros(1,length(t));
    W1initial = 0;
    W2initial = 0;
    W1d1 = zeros(1,length(t));    
    W2d1 = zeros(1,length(t));
    dW1d1 = zeros(1,length(t));
    dW2d1 = zeros(1,length(t));
    W1d2 = zeros(1,length(t));    
    W2d2 = zeros(1,length(t));
    dW1d2 = zeros(1,length(t));
    dW2d2 = zeros(1,length(t));

    % Form the increments:
    for n = 1:length(t)-1
        dW1(n) = sqrt(t(n+1)-t(n))*randn(1,1);
        dW2(n) = sqrt(t(n+1)-t(n))*randn(1,1);
    end

    % Form the processes:
    for n = 1:length(t)
        if n == 1
            W1(1,n) = dW1initial(1,1);
            W2(1,n) = dW2initial(1,1);
        else
            W1(1,n) = W1(1,n-1) + dW1(1,n-1);
            W2(1,n) = W2(1,n-1) + dW2(1,n-1);
        end
    end

    % Form Wjdk(t)=Wj(t-dk):   
    for n = 1:length(t)
        if t(n) == t(td1index(1))
            W1d1(1,n) = W1initial;
            W2d1(1,n) = W2initial;
        elseif t(n) > d1
            j = find(t==t(n)-d1);
            W1d1(1,n) = W1(1,j);
            W2d1(1,n) = W2(1,j);
        end
        if t(n) == t(td2index(1))
            W1d2(1,n) = W1initial;
            W2d2(1,n) = W2initial;
        elseif t(n) > d2
            j = find(t==t(n)-d2);
            W1d2(1,n) = W1(1,j);
            W2d2(1,n) = W2(1,j);
        end
    end
    for n = td1index(1):length(t)-1
        dW1d1(1,n) = W1d1(1,n+1)-W1d1(1,n);
        dW2d1(1,n) = W2d1(1,n+1)-W2d1(1,n);
    end
    for n = td2index(1):length(t)-1
        dW1d2(1,n) = W1d2(1,n+1)-W1d2(1,n);
        dW2d2(1,n) = W2d2(1,n+1)-W2d2(1,n);
    end
    
    % To test these, use the following lines in the command window:
    % plot([0,t],[0,W1],'k',[0,t],[0,W1d1],'b');legend('$W_1$','$W_1^{\tau_1}$','Interpreter','latex');if max(W1(1,1:find(t==T-d1)) == W1d1(td1index(1)+1:td1index(1)+length(W1(1,1:find(t==T-d1))))) == 1 && min(W1(1,1:find(t==T-d1)) == W1d1(td1index(1)+1:td1index(1)+length(W1(1,1:find(t==T-d1))))) == 1; disp('W1 and W1d1 calculated correctly.'); else disp('W1 and W1d1 calculated incorrectly.'); end; 
    % plot([0,t],[0,W1],'k',[0,t],[0,W1d2],'b');legend('$W_1$','$W_1^{\tau_2}$','Interpreter','latex');if max(W1(1,1:find(t==T-d2)) == W1d2(td2index(1)+1:td2index(1)+length(W1(1,1:find(t==T-d2))))) == 1 && min(W1(1,1:find(t==T-d2)) == W1d2(td2index(1)+1:td2index(1)+length(W1(1,1:find(t==T-d2))))) == 1; disp('W1 and W1d2 calculated correctly.'); else disp('W1 and W1d2 calculated incorrectly.'); end; 
    % plot([0,t],[0,W2],'k',[0,t],[0,W2d1],'b');legend('$W_2$','$W_2^{\tau_1}$','Interpreter','latex');if max(W2(1,1:find(t==T-d1)) == W2d1(td1index(1)+1:td1index(1)+length(W2(1,1:find(t==T-d1))))) == 1 && min(W2(1,1:find(t==T-d1)) == W2d1(td1index(1)+1:td1index(1)+length(W2(1,1:find(t==T-d1))))) == 1; disp('W2 and W2d1 calculated correctly.'); else disp('W2 and W2d1 calculated incorrectly.'); end; 
    % plot([0,t],[0,W2],'k',[0,t],[0,W2d2],'b');legend('$W_2$','$W_2^{\tau_2}$','Interpreter','latex');if max(W2(1,1:find(t==T-d2)) == W2d2(td2index(1)+1:td2index(1)+length(W2(1,1:find(t==T-d2))))) == 1 && min(W2(1,1:find(t==T-d2)) == W2d2(td2index(1)+1:td2index(1)+length(W2(1,1:find(t==T-d2))))) == 1; disp('W2 and W2d2 calculated correctly.'); else disp('W2 and W2d2 calculated incorrectly.'); end; 

%% SIMULATION OF REFERENCE SOLUTION:

    % Construct the reference (analytic) solution:
    Xref = zeros(2,length(t));
    % Initial step:
    Xn = X0; Xd1 = phi(-d1); Xd2 = phi(-d2);
    a = A0*Xn + B0(Xn,Xd1,Xd2);
    b1 = A1*Xn + B1(Xn,Xd1,Xd2);  b1x = A1+B1x(Xn,Xd1,Xd2);
    b2 = A2*Xn + B2(Xn,Xd1,Xd2);  b2x = A2+B2x(Xn,Xd1,Xd2);
    I11 = (dW1initial(1,1)^2-t(1))/2; I22 = (dW2initial(1,1)^2-t(1))/2;
        % I12ref calculation:
        I12 = dW1initial(1,1)*dW2initial(1,1)/2;
        I21 = dW1initial(1,1)*dW2initial(1,1) - I12;
    % Update these Iijs for these refined ones.
    Lx1b1 = (b1x)*b1;   Lx2b2 = (b2x)*b2;   Lx1b2 = (b2x)*b1;   Lx2b1 = (b1x)*b2;
    Xref(:,1) = Xn + a*t(1) + b1*dW1initial(1,1) + b2*dW2initial(1,1) + Lx1b1*I11 + Lx2b1*I21 + Lx1b2*I12 + Lx2b2*I22;
    % Inductive steps:
    for n = 1:length(t)-1
        Xn = Xref(:,n);
        if t(n) <= d1 && t(n) <= d2
            Xd1 = phi(t(n)-d1); Xd1d1 = phi(t(n)-d1-d1); Xd2 = phi(t(n)-d2); Xd2d2 = phi(t(n)-d2-d2);
        elseif d1 < t(n) && t(n) <= 2*d1 && t(n) <= d2
            i1 = find(t==t(n)-d1);
            Xd1 = Xref(:,i1); Xd1d1 = phi(t(n)-d1-d1); Xd2 = phi(t(n)-d2); Xd2d2 = phi(t(n)-d2-d2);
        elseif t(n) > 2*d1 && t(n) <= d2
            i1 = find(t==t(n)-d1); i2 = find(t==t(n)-2*d1);
            Xd1 = Xref(:,i1); Xd1d1 = Xref(:,i2); Xd2 = phi(t(n)-d2); Xd2d2 = phi(t(n)-2*d2);
        elseif t(n) <= d1 && d2 < t(n) && t(n) <= 2*d2
            j1 = find(t==t(n)-d2);
            Xd1 = phi(t(n)-d1); Xd1d1 = phi(t(n)-2*d1); Xd2 = Xref(:,j1); Xd2d2 = phi(t(n)-2*d2);
        elseif d1 < t(n) && t(n) <= 2*d1 && d2 < t(n) && t(n) <= 2*d2
            i1 = find(t==t(n)-d1); j1 = find(t==t(n)-d2);
            Xd1 = Xref(:,i1); Xd1d1 = phi(t(n)-2*d1); Xd2 = Xref(:,j1); Xd2d2 = phi(t(n)-2*d2);
        elseif t(n) > 2*d1 && d2 < t(n) && t(n) <= 2*d2
            i1 = find(t==t(n)-d1); i2 = find(t==t(n)-2*d1); j1 = find(t==t(n)-d2);
            Xd1 = Xref(:,i1); Xd1d1 = Xref(:,i2); Xd2 = Xref(:,j1); Xd2d2 = phi(t(n)-2*d2);
        elseif t(n) <= d1 && 2*d2 < t(n)
            j1 = find(t==t(n)-d2); j2 = find(t==t(n)-2*d2);
            Xd1 = phi(t(n)-d1); Xd1d1 = phi(t(n)-d1-d1); Xd2 = Xref(:,j1); Xd2d2 = Xref(:,j2);
        elseif d1 < t(n) && t(n) <= 2*d1 && 2*d2 < t(n)
            i1 = find(t==t(n)-d1); j1 = find(t==t(n)-d2); j2 = find(t==t(n)-2*d2);
            Xd1 = Xref(:,i1); Xd1d1 = phi(t(n)-d1-d1); Xd2 = Xref(:,j1); Xd2d2 = Xref(:,j2);            
        elseif t(n) > 2*d1 && 2*d2 < t(n)
            i1 = find(t==t(n)-d1); i2 = find(t==t(n)-2*d1); j1 = find(t==t(n)-d2); j2 = find(t==t(n)-2*d2);
            Xd1 = Xref(:,i1); Xd1d1 = Xref(:,i2); Xd2 = Xref(:,j1); Xd2d2 = Xref(:,j2);
        end
        if t(n) <= d1 + d2
            Xd1d2 = phi(t(n)-d1-d2);
        elseif t(n) > d1 + d2
            i1j1 = find(t==t(n)-d1-d2); Xd1d2 = Xref(:,i1j1);
        end
        Xd2d1 = Xd1d2;
        hsize = t(n+1)-t(n);
        a = A0*Xn + B0(Xn,Xd1,Xd2);
        b1 = A1*Xn + B1(Xn,Xd1,Xd2);  b1x = A1+B1x(Xn,Xd1,Xd2);
        b2 = A2*Xn + B2(Xn,Xd1,Xd2);  b2x = A2+B2x(Xn,Xd1,Xd2);
        I11Ref = (dW1(1,n)^2-hsize)/2; I22Ref = (dW2(1,n)^2-hsize)/2; I12Ref = dW1(1,n)*dW2(1,n)/2; I21Ref = dW1(1,n)*dW2(1,n) - I12Ref;
        Lx1b1 = (b1x)*b1;   Lx2b2 = (b2x)*b2;   Lx1b2 = (b2x)*b1;   Lx2b1 = (b1x)*b2;
        b1xd1 = B1xd1(Xn,Xd1,Xd2); b2xd1 = B2xd1(Xn,Xd1,Xd2); b1d1 = A1*Xd1+ B1(Xd1,Xd1d1,Xd2d1); b2d1 = A2*Xd1+ B2(Xd1,Xd1d1,Xd2d1);
        Lxd_one1b1 = (b1xd1)*b1d1; Lxd_one2b2 = (b2xd1)*b2d1;   Lxd_one1b2 = (b2xd1)*b1d1;   Lxd_one2b1 = (b1xd1)*b2d1;
        if t(n) < d1
            I11d1Ref = 0; I12d1Ref = 0; I21d1Ref = 0; I22d1Ref = 0;
        else
            I11d1Ref = dW1d1(1,n)*dW1(1,n)/2; I12d1Ref = dW1d1(1,n)*dW2(1,n)/2; I21d1Ref = dW1(1,n)*dW2d1(1,n)/2; I22d1Ref = dW2(1,n)*dW2d1(1,n)/2;
        end

        b1xd2 = B1xd2(Xn,Xd1,Xd2); b2xd2 = B2xd2(Xn,Xd1,Xd2); b1d2 = A1*Xd2+ B1(Xd2,Xd1d2,Xd2d2); b2d2 = A2*Xd2+ B2(Xd2,Xd1d2,Xd2d2);
        Lxd_two1b1 = (b1xd2)*b1d2; Lxd_two2b2 = (b2xd2)*b2d2;   Lxd_two1b2 = (b2xd2)*b1d2;   Lxd_two2b1 = (b1xd2)*b2d2;
        if t(n) < d2
            I11d2Ref = 0; I12d2Ref = 0; I21d2Ref = 0; I22d2Ref = 0;
        else
            I11d2Ref = dW1d2(1,n)*dW1(1,n)/2; I12d2Ref = dW1d2(1,n)*dW2(1,n)/2; I21d2Ref = dW1(1,n)*dW2d2(1,n)/2; I22d2Ref = dW2(1,n)*dW2d2(1,n)/2;
        end

        Xref(:,n+1) = Xn + a*hsize + b1*dW1(1,n) + b2*dW2(1,n) + Lx1b1*I11Ref + Lx2b1*I21Ref + Lx1b2*I12Ref + Lx2b2*I22Ref + Lxd_one1b1*I11d1Ref + Lxd_one1b2*I12d1Ref + Lxd_one2b1*I21d1Ref + Lxd_one2b2*I22d1Ref  +  Lxd_two1b1*I11d2Ref + Lxd_two1b2*I12d2Ref + Lxd_two2b1*I21d2Ref + Lxd_two2b2*I22d2Ref;
    end

    for ti = 1:NO
        Xrefvalues(:,ti,trial) = Xref(:,tindexobservations(ti));
    end
    if trial == samplepath
        sampleXref = Xref(samplecomponent,:);
    end

%% SCHEMES:

    for h = 1:length(hvalues)

        % Define some constants that are used throughout the following part:
        step = hvalues(h)/refdt;    % This is the index multiplier, for how many index positions comprise a single step of the simulations.  The idea is that t_n = t(n*step).
        Delt = hvalues(h);          % Step size.
        hsize = Delt;   % Just another name for it.
        A0hatDelt = (A0 - A1^2/2 - A2^2/2)*Delt;
        schememesh = mesh0(step):mesh0(step):T;
        % TOscheme = [schememesh,(1:N1)*d1,(1:N2)*d2,tobservations];
        % TCscheme = TOscheme;
        % for i = 0:N1
        %     for j = 0:N2
        %         TCscheme = [TCscheme,TCscheme - i*d1 - j*d2];
        %         TCscheme = sort(TCscheme);    TCscheme = TCscheme(0<TCscheme&TCscheme<=T);    TCscheme = unique(TCscheme);
        %     end
        % end
        % TCscheme = sort(TCscheme);    TCscheme = TCscheme(0<TCscheme&TCscheme<=T);    TCscheme = unique(TCscheme);
        % 
        % schememesh = TCscheme; % Time for the schemes.
        N = length(schememesh);           % This is the number of steps that are made in [0,T].

        % Find the index positions in schememesh to observe the error graphs at:
        schememeshindexobservations = zeros(1,length(tobservations));
        for i = 1:length(tobservations)
            schememeshindexobservations(i) = find(schememesh==tobservations(i));
        end

        % Initialise DelW1, DelW2, DelW1d, and DelW2d:
        DelW1 = zeros(1,N-1);
        DelW2 = zeros(1,N-1);
        DelW1d1 = zeros(1,N-1);
        DelW2d1 = zeros(1,N-1);
        DelW1d2 = zeros(1,N-1);
        DelW2d2 = zeros(1,N-1);
        i = find(t==schememesh(1));
        DelW1initial = W1(1,i)-W1initial;
        DelW2initial = W2(1,i)-W2initial;
        DelW1d1initial = zeros(1,1);
        DelW2d1initial = zeros(1,1);
        DelW1d2initial = zeros(1,1);
        DelW2d2initial = zeros(1,1);

        % Form DelW1 and DelW2:
        for n = 1:N-1
            i1 = find(t==schememesh(n+1));
            i2 = find(t==schememesh(n));
            DelW1(1,n) = W1(1,i1)-W1(1,i2);
            DelW2(1,n) = W2(1,i1)-W2(1,i2);
        end
        % Form DelW1d1 and DelW2d1:
        for n = 1:N-1
            if schememesh(n) == d1
                j1 = find(t==schememesh(n+1)-d1);
                j2 = find(t==schememesh(n)-d1);
                DelW1d1(1,n) = W1(1,j1) - W1initial(1,1);
                DelW2d1(1,n) = W2(1,j1) - W2initial(1,1);
            elseif schememesh(n) > d1
                j1 = find(t==schememesh(n+1)-d1);
                j2 = find(t==schememesh(n)-d1);
                DelW1d1(1,n) = W1(1,j1) - W1(1,j2);
                DelW2d1(1,n) = W2(1,j1) - W2(1,j2);
            end
        end
        % Form DelW1d2 and DelW2d2:
        for n = 1:N-1
            if schememesh(n) == d2
                j1 = find(t==schememesh(n+1)-d2);
                j2 = find(t==schememesh(n)-d2);
                DelW1d2(1,n) = W1(1,j1) - W1initial(1,1);
                DelW2d2(1,n) = W2(1,j1) - W2initial(1,1);
            elseif schememesh(n) > d2
                j1 = find(t==schememesh(n+1)-d2);
                j2 = find(t==schememesh(n)-d2);
                DelW1d2(1,n) = W1(1,j1) - W1(1,j2);
                DelW2d2(1,n) = W2(1,j1) - W2(1,j2);
            end
        end

        % To test this in the command window:
        % hold on; plot([0,t],[0,W1],'k'); plot([0,t],[0,W1d],'b'); plot([schememesh(1:end-1)],[0,cumsum(DelW1d(1:end-1))],'r:'); legend('$W_1$','$W_1^{\tau}$','$W_1^{\tau}$ from $\Delta W_1^{\tau}$','Interpreter','latex')


        % Initialise the schemes:
        EM = zeros(2,N);       % Euler--Maruyama 
        MilSim = zeros(2,N);   % Milstein (Simple)
        MEM = zeros(2,N);      % Magnus--EM
        MMSim = zeros(2,N);    % Magnus--Mistein (Simple)        
        MilRef = zeros(2,N);   % Milstein (Refined)
        MMRef = zeros(2,N);    % Magnus--Milstein (Refined)
        
        % The first step, from t=0 to t=h, using initial increments, such
        % as DelW1initial.  Effectively, this is n = 0:

    % INITIAL STEP:

        % EM:
        Yn = X0; Yd1 = phi(-d1); Yd2 = phi(-d2);
        a = A0*Yn + B0(Yn,Yd1,Yd2);
        b1 = A1*Yn + B1(Yn,Yd1,Yd2);
        b2 = A2*Yn + B2(Yn,Yd1,Yd2);
        EM(:,1) = X0 + a*schememesh(1) + b1*DelW1initial(1,1) + b2*DelW2initial(1,1);

        % MilSim:
        Yn = X0; Yd1 = phi(-d1); Yd2 = phi(-d2);
        a = A0*Yn + B0(Yn,Yd1,Yd2);
        b1 = A1*Yn + B1(Yn,Yd1,Yd2); b1x = A1+B1x(Yn,Yd1,Yd2);
        b2 = A2*Yn + B2(Yn,Yd1,Yd2); b2x = A2+B2x(Yn,Yd1,Yd2);
        I11Sim = (DelW1initial(1,1)^2-schememesh(1))/2; I22Sim = (DelW2initial(1,1)^2-schememesh(1))/2; I12Sim = DelW1initial(1,1)*DelW2initial(1,1)/2; I21Sim = DelW1initial(1,1)*DelW2initial(1,1)/2;
        Lx1b1 = (b1x)*b1;   Lx2b2 = (b2x)*b2;   Lx1b2 = (b2x)*b1;   Lx2b1 = (b1x)*b2;
        MilSim(:,1) = Yn + a*schememesh(1) + b1*DelW1initial(1,1) + b2*DelW2initial(1,1) + Lx1b1*I11Sim + Lx2b1*I21Sim + Lx1b2*I12Sim + Lx2b2*I22Sim;

        % MEM:
        Yn = X0; Yd1 = phi(-d1); Yd2 = phi(-d2);
        a = B0(Yn,Yd1,Yd2);
        b1 = B1(Yn,Yd1,Yd2);
        b2 = B2(Yn,Yd1,Yd2);
        expOmega1 = expm( A0hat*schememesh(1) + A1*DelW1initial(1,1) + A2*DelW2initial(1,1) );
        MEM(:,1) = expOmega1*(Yn + (a - A1*b1 - A2*b2)*schememesh(1) + b1*DelW1initial(1,1) + b2*DelW2initial(1,1) );

        % MMSim:
        Yn = X0; Yd1 = phi(-d1); Yd2 = phi(-d2);
        a = B0(Yn,Yd1,Yd2);
        b1 = B1(Yn,Yd1,Yd2); b1x = B1x(Yn,Yd1,Yd2);
        b2 = B2(Yn,Yd1,Yd2); b2x = B2x(Yn,Yd1,Yd2);
        I10Sim = DelW1initial(1,1)*schememesh(1)/2; I01Sim = DelW1initial(1,1)*schememesh(1) - I10Sim; I20Sim = DelW2initial(1,1)*schememesh(1)/2; I02Sim = DelW2initial(1,1)*schememesh(1) - I20Sim;
        expOmega2 = expm( A0hat*schememesh(1) + A1*DelW1initial(1,1) + A2*DelW2initial(1,1) + 1/2*( LieA0A1*(I10Sim-I01Sim)+LieA0A2*(I20Sim-I02Sim)+LieA1A2*(I21Sim-I12Sim) ) );
        Hx1b1 = (b1x)*(A1*Yn+b1)-A1*b1; % This 0 is actually the Jacobian, d(b1)/dx.   j = 1, l = 1
        Hx1b2 = (b1x)*(A2*Yn+b2)-A2*b1; %  j = 1, l = 2
        Hx2b1 = (b2x)*(A1*Yn+b1)-A1*b2; %  j = 2, l = 1
        Hx2b2 = (b2x)*(A2*Yn+b2)-A2*b2; %  j = 2, l = 1
        MMSim(:,1) = expOmega2*(Yn + (a - A1*b1 - A2*b2)*schememesh(1) + b1*DelW1initial(1,1) + b2*DelW2initial(1,1) + Hx1b1*I11Sim + Hx1b2*I21Sim + Hx2b1*I12Sim + Hx2b2*I22Sim );

        % Refined schemes:

            % MilRef:
            Yn = X0; Yd1 = phi(-d1); Yd2 = phi(-d2);
            a = A0*Yn + B0(Yn,Yd1,Yd2);
            b1 = A1*Yn + B1(Yn,Yd1,Yd2); b1x = A1+B1x(Yn,Yd1,Yd2);
            b2 = A2*Yn + B2(Yn,Yd1,Yd2); b2x = A2+B2x(Yn,Yd1,Yd2);
            I11Ref = (DelW1initial(1,1)^2-schememesh(1))/2; I22Ref = (DelW2initial(1,1)^2-schememesh(1))/2;
                % I12ref calculation:
                I12Ref = dW1initial(1,1)*dW2initial(1,1)/2;
                for j=1:find(t==schememesh(1))-1
                    I12Ref = I12Ref + dW1(1,j)*dW2(1,j)/2 + dW2(1,j)*(W1(1,j)-W1initial(1,1));
                end
                I21Ref = DelW1initial(1,1)*DelW2initial(1,1) - I12Ref;
            Lx1b1 = (b1x)*b1;   Lx2b2 = (b2x)*b2;   Lx1b2 = (b2x)*b1;   Lx2b1 = (b1x)*b2;
            MilRef(:,1) = Yn + a*schememesh(1) + b1*DelW1initial(1,1) + b2*DelW2initial(1,1) + Lx1b1*I11Ref + Lx2b1*I21Ref + Lx1b2*I12Ref + Lx2b2*I22Ref;

            % MMRef:
            Yn = X0; Yd1 = phi(-d1); Yd2 = phi(-d2);
            a = B0(Yn,Yd1,Yd2);
            b1 = B1(Yn,Yd1,Yd2); b1x = B1x(Yn,Yd1,Yd2);
            b2 = B2(Yn,Yd1,Yd2); b2x = B2x(Yn,Yd1,Yd2);
            I10Ref = dW1initial(1,1)*(t(1)-0)/2;
                for j=1:find(t==schememesh(1))-1
                    I10Ref = I10Ref + dW1(1,j)*(t(1)-0)/2 + refdt*(W1(1,j)-W1initial(1,1));
                end
            I01Ref = DelW1initial(1,1)*schememesh(1) - I10Ref;
            I20Ref = dW2initial(1,1)*(t(1)-0)/2;
                for j=1:find(t==schememesh(1))-1
                    I20Ref = I20Ref + dW2(1,j)*(t(1)-0)/2 + refdt*(W2(1,j)-W2initial(1,1));
                end
            I02Ref = DelW2initial(1,1)*schememesh(1) - I20Ref;
            expOmega2 = expm( A0hat*schememesh(1) + A1*DelW1initial(1,1) + A2*DelW2initial(1,1) + 1/2*( LieA0A1*(I10Ref-I01Ref)+LieA0A2*(I20Ref-I02Ref)+LieA1A2*(I21Ref-I12Ref) ) );
            Hx1b1 = (b1x)*(A1*Yn+b1)-A1*b1; % This 0 is actually the Jacobian, d(b1)/dx.   j = 1, l = 1
            Hx1b2 = (b1x)*(A2*Yn+b2)-A2*b1; %  j = 1, l = 2
            Hx2b1 = (b2x)*(A1*Yn+b1)-A1*b2; %  j = 2, l = 1
            Hx2b2 = (b2x)*(A2*Yn+b2)-A2*b2; %  j = 2, l = 1
            MMRef(:,1) = expOmega2*(Yn + (a - A1*b1 - A2*b2)*schememesh(1) + b1*DelW1initial(1,1) + b2*DelW2initial(1,1) + Hx1b1*I11Ref + Hx1b2*I21Ref + Hx2b1*I12Ref + Hx2b2*I22Ref );

    % ITERATIVE STEPS:

        for n = 1:N-1

            tn = schememesh(n);

            % Set the schememesh indices:
            i1f = find(t==schememesh(n)); % f for "finest" mesh.
            i2f = find(t==schememesh(n+1)); % f for "finest" mesh.AAAA

            % Set the delay indices:
            tnd1pre_index = find(schememesh <= schememesh(n)-d1,1,'last'); tnd1pre = schememesh(tnd1pre_index); if isempty(tnd1pre); tnd1pre = 0; end
                tnd1post_index = find(schememesh > schememesh(n)-d1,1,'first'); tnd1post = schememesh(tnd1post_index); if isempty(tnd1post); tnd1post = 0; end
            tnd2pre_index = find(schememesh <= schememesh(n)-d2,1,'last'); tnd2pre = schememesh(tnd2pre_index); if isempty(tnd2pre); tnd2pre = 0; end
                tnd2post_index = find(schememesh > schememesh(n)-d2,1,'first'); tnd2post = schememesh(tnd2post_index); if isempty(tnd2post); tnd2post = 0; end
            tnd1d1pre_index = find(schememesh <= schememesh(n)-d1-d1,1,'last'); tnd1d1pre = schememesh(tnd1d1pre_index); if isempty(tnd1d1pre); tnd1d1pre = 0; end
                tnd1d1post_index = find(schememesh > schememesh(n)-d1-d1,1,'first'); tnd1d1post = schememesh(tnd1d1post_index); if isempty(tnd1d1post); tnd1d1post = 0; end
            tnd1d2pre_index = find(schememesh <= schememesh(n)-d1-d2,1,'last'); tnd1d2pre = schememesh(tnd1d2pre_index); if isempty(tnd1d2pre); tnd1d2pre = 0; end
                tnd1d2post_index = find(schememesh > schememesh(n)-d1-d2,1,'first'); tnd1d2post = schememesh(tnd1d2post_index); if isempty(tnd1d2post); tnd1d2post = 0; end
            tnd2d2pre_index = find(schememesh <= schememesh(n)-d2-d2,1,'last'); tnd2d2pre = schememesh(tnd2d2pre_index); if isempty(tnd2d2pre); tnd2d2pre = 0; end
                tnd2d2post_index = find(schememesh > schememesh(n)-d2-d2,1,'first'); tnd2d2post = schememesh(tnd2d2post_index); if isempty(tnd2d2post); tnd2d2post = 0; end
        
            % EM:
            Yn = EM(:,n);
            if schememesh(n) <= d1 && schememesh(n) <= d2
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = EM(:,tnd1pre_index); end; Yd1post = EM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = EM(:,tnd1pre_index); end; Yd1post = EM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = EM(:,tnd1d1pre_index); end; Yd1d1post = EM(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = EM(:,tnd2pre_index); end; Yd2post = EM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = EM(:,tnd1pre_index); end; Yd1post = EM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = EM(:,tnd2pre_index); end; Yd2post = EM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = EM(:,tnd1pre_index); end; Yd1post = EM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = EM(:,tnd1d1pre_index); end; Yd1d1post = EM(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = EM(:,tnd2pre_index); end; Yd2post = EM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && 2*d2 < schememesh(n)
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = EM(:,tnd2pre_index); end; Yd2post = EM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = EM(:,tnd2d2pre_index); end; Yd2d2post = EM(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = EM(:,tnd1pre_index); end; Yd1post = EM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = EM(:,tnd2pre_index); end; Yd2post = EM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = EM(:,tnd2d2pre_index); end; Yd2d2post = EM(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif schememesh(n) > 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = EM(:,tnd1pre_index); end; Yd1post = EM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = EM(:,tnd1d1pre_index); end; Yd1d1post = EM(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = EM(:,tnd2pre_index); end; Yd2post = EM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = EM(:,tnd2d2pre_index); end; Yd2d2post = EM(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
            end
            if schememesh(n) <= d1 + d2
                Yd1d2 = phi(schememesh(n)-d1-d2);
            elseif schememesh(n) > d1 + d2
                if isempty(tnd1d2pre_index); Yd1d2pre = phi(0); else; Yd1d2pre = EM(:,tnd1d2pre_index); end; Yd1d2post = EM(:,tnd1d2post_index);
                    Yd1d2 = Yd1d2pre + (Yd1d2post-Yd1d2pre)*((tn-d1-d2)-tnd1d2pre)/hsize;
            end
            Yd2d1 = Yd1d2;
            a = A0*Yn + B0(Yn,Yd1,Yd2);
            b1 = A1*Yn + B1(Yn,Yd1,Yd2);
            b2 = A2*Yn + B2(Yn,Yd1,Yd2);
            EM(:,n+1) = Yn + a*hsize + b1*DelW1(1,n) + b2*DelW2(1,n);

            % MilSim:
            Yn = MilSim(:,n);
            if schememesh(n) <= d1 && schememesh(n) <= d2
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilSim(:,tnd1pre_index); end; Yd1post = MilSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilSim(:,tnd1pre_index); end; Yd1post = MilSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MilSim(:,tnd1d1pre_index); end; Yd1d1post = MilSim(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilSim(:,tnd2pre_index); end; Yd2post = MilSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilSim(:,tnd1pre_index); end; Yd1post = MilSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilSim(:,tnd2pre_index); end; Yd2post = MilSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilSim(:,tnd1pre_index); end; Yd1post = MilSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MilSim(:,tnd1d1pre_index); end; Yd1d1post = MilSim(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilSim(:,tnd2pre_index); end; Yd2post = MilSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && 2*d2 < schememesh(n)
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilSim(:,tnd2pre_index); end; Yd2post = MilSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MilSim(:,tnd2d2pre_index); end; Yd2d2post = MilSim(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilSim(:,tnd1pre_index); end; Yd1post = MilSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilSim(:,tnd2pre_index); end; Yd2post = MilSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MilSim(:,tnd2d2pre_index); end; Yd2d2post = MilSim(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif schememesh(n) > 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilSim(:,tnd1pre_index); end; Yd1post = MilSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MilSim(:,tnd1d1pre_index); end; Yd1d1post = MilSim(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilSim(:,tnd2pre_index); end; Yd2post = MilSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MilSim(:,tnd2d2pre_index); end; Yd2d2post = MilSim(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
            end
            if schememesh(n) <= d1 + d2
                Yd1d2 = phi(schememesh(n)-d1-d2);
            elseif schememesh(n) > d1 + d2
                if isempty(tnd1d2pre_index); Yd1d2pre = phi(0); else; Yd1d2pre = MilSim(:,tnd1d2pre_index); end; Yd1d2post = MilSim(:,tnd1d2post_index);
                    Yd1d2 = Yd1d2pre + (Yd1d2post-Yd1d2pre)*((tn-d1-d2)-tnd1d2pre)/hsize;
            end
            Yd2d1 = Yd1d2;
            %a,b1,b2,...
            a = A0*Yn + B0(Yn,Yd1,Yd2);
            b1 = A1*Yn + B1(Yn,Yd1,Yd2); b1x = A1+B1x(Yn,Yd1,Yd2);
            b2 = A2*Yn + B2(Yn,Yd1,Yd2); b2x = A2+B2x(Yn,Yd1,Yd2);
            I11Sim = (DelW1(1,n)^2-hsize)/2; I22Sim = (DelW2(1,n)^2-hsize)/2; I12Sim = DelW1(1,n)*DelW2(1,n)/2; I21Sim = DelW1(1,n)*DelW2(1,n)/2;
            Lx1b1 = (b1x)*b1;   Lx2b2 = (b2x)*b2;   Lx1b2 = (b2x)*b1;   Lx2b1 = (b1x)*b2;
            b1xd1 = B1xd1(Yn,Yd1,Yd2); b2xd1 = B2xd1(Yn,Yd1,Yd2); b1d1 = A1*Yd1 + B1(Yd1,Yd1d1,Yd2d1); b2d1 = A2*Yd1 + B2(Yd1,Yd1d1,Yd2d1);
            Lxd_one1b1 = (b1xd1)*b1d1; Lxd_one2b2 = (b2xd1)*b2d1;   Lxd_one1b2 = (b2xd1)*b1d1;   Lxd_one2b1 = (b1xd1)*b2d1;
            if schememesh(n) < d1
                I11d1Sim = 0; I22d1Sim = 0; I12d1Sim = 0; I21d1Sim = 0;
            else
                I11d1Sim = DelW1d1(1,n)*DelW1(1,n)/2; I22d1Sim = DelW2d1(1,n)*DelW2(1,n)/2; I12d1Sim = DelW1d1(1,n)*DelW2(1,n)/2; I21d1Sim = DelW2d1(1,n)*DelW1(1,n)/2;
            end
            b1xd2 = B1xd2(Yn,Yd1,Yd2); b2xd2 = B2xd2(Yn,Yd1,Yd2); b1d2 = A1*Yd2 + B1(Yd2,Yd1d2,Yd2d2); b2d2 = A2*Yd2 + B2(Yd2,Yd1d2,Yd2d2);
            Lxd_two1b1 = (b1xd2)*b1d2; Lxd_two2b2 = (b2xd2)*b2d2;   Lxd_two1b2 = (b2xd2)*b1d2;   Lxd_two2b1 = (b1xd2)*b2d2;
            if schememesh(n) < d2
                I11d2Sim = 0; I22d2Sim = 0; I12d2Sim = 0; I21d2Sim = 0;
            else
                I11d2Sim = DelW1d2(1,n)*DelW1(1,n)/2; I22d2Sim = DelW2d2(1,n)*DelW2(1,n)/2; I12d2Sim = DelW1d2(1,n)*DelW2(1,n)/2; I21d2Sim = DelW2d2(1,n)*DelW1(1,n)/2;
            end
            MilSim(:,n+1) = Yn + a*hsize + b1*DelW1(1,n) + b2*DelW2(1,n) + Lx1b1*I11Sim + Lx1b2*I12Sim + Lx2b1*I21Sim + Lx2b2*I22Sim + Lxd_one1b1*I11d1Sim + Lxd_one1b2*I12d1Sim + Lxd_one2b1*I21d1Sim + Lxd_one2b2*I22d1Sim  +  Lxd_two1b1*I11d2Sim + Lxd_two1b2*I12d2Sim + Lxd_two2b1*I21d2Sim + Lxd_two2b2*I22d2Sim;

            % MEM:
            Yn = MEM(:,n);
            if schememesh(n) <= d1 && schememesh(n) <= d2
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MEM(:,tnd1pre_index); end; Yd1post = MEM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MEM(:,tnd1pre_index); end; Yd1post = MEM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MEM(:,tnd1d1pre_index); end; Yd1d1post = MEM(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MEM(:,tnd2pre_index); end; Yd2post = MEM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MEM(:,tnd1pre_index); end; Yd1post = MEM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MEM(:,tnd2pre_index); end; Yd2post = MEM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MEM(:,tnd1pre_index); end; Yd1post = MEM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MEM(:,tnd1d1pre_index); end; Yd1d1post = MEM(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MEM(:,tnd2pre_index); end; Yd2post = MEM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && 2*d2 < schememesh(n)
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MEM(:,tnd2pre_index); end; Yd2post = MEM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MEM(:,tnd2d2pre_index); end; Yd2d2post = MEM(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MEM(:,tnd1pre_index); end; Yd1post = MEM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MEM(:,tnd2pre_index); end; Yd2post = MEM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MEM(:,tnd2d2pre_index); end; Yd2d2post = MEM(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif schememesh(n) > 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MEM(:,tnd1pre_index); end; Yd1post = MEM(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MEM(:,tnd1d1pre_index); end; Yd1d1post = MEM(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MEM(:,tnd2pre_index); end; Yd2post = MEM(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MEM(:,tnd2d2pre_index); end; Yd2d2post = MEM(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
            end
            if schememesh(n) <= d1 + d2
                Yd1d2 = phi(schememesh(n)-d1-d2);
            elseif schememesh(n) > d1 + d2
                if isempty(tnd1d2pre_index); Yd1d2pre = phi(0); else; Yd1d2pre = MEM(:,tnd1d2pre_index); end; Yd1d2post = MEM(:,tnd1d2post_index);
                    Yd1d2 = Yd1d2pre + (Yd1d2post-Yd1d2pre)*((tn-d1-d2)-tnd1d2pre)/hsize;
            end
            Yd2d1 = Yd1d2;
            a = B0(Yn,Yd1,Yd2);
            b1 = B1(Yn,Yd1,Yd2);
            b2 = B2(Yn,Yd1,Yd2);
            expOmega1 = expm( A0hat*hsize + A1*DelW1(1,n) + A2*DelW2(1,n) );
            MEM(:,n+1) = expOmega1*(Yn + (a - A1*b1 - A2*b2)*hsize + b1*DelW1(1,n) + b2*DelW2(1,n) );

            % MMSim:
            Yn = MMSim(:,n);
            if schememesh(n) <= d1 && schememesh(n) <= d2
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMSim(:,tnd1pre_index); end; Yd1post = MMSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && schememesh(n) <= d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMSim(:,tnd1pre_index); end; Yd1post = MMSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MMSim(:,tnd1d1pre_index); end; Yd1d1post = MMSim(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMSim(:,tnd2pre_index); end; Yd2post = MMSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMSim(:,tnd1pre_index); end; Yd1post = MMSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMSim(:,tnd2pre_index); end; Yd2post = MMSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) > 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMSim(:,tnd1pre_index); end; Yd1post = MMSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MMSim(:,tnd1d1pre_index); end; Yd1d1post = MMSim(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMSim(:,tnd2pre_index); end; Yd2post = MMSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                Yd2d2 = phi(schememesh(n)-d2-d2);
            elseif schememesh(n) <= d1 && 2*d2 < schememesh(n)
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMSim(:,tnd2pre_index); end; Yd2post = MMSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MMSim(:,tnd2d2pre_index); end; Yd2d2post = MMSim(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMSim(:,tnd1pre_index); end; Yd1post = MMSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMSim(:,tnd2pre_index); end; Yd2post = MMSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MMSim(:,tnd2d2pre_index); end; Yd2d2post = MMSim(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                Yd1d1 = phi(schememesh(n)-d1-d1);
            elseif schememesh(n) > 2*d1 && 2*d2 < schememesh(n)
                if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMSim(:,tnd1pre_index); end; Yd1post = MMSim(:,tnd1post_index);
                    Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MMSim(:,tnd1d1pre_index); end; Yd1d1post = MMSim(:,tnd1d1post_index);
                    Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMSim(:,tnd2pre_index); end; Yd2post = MMSim(:,tnd2post_index);
                    Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MMSim(:,tnd2d2pre_index); end; Yd2d2post = MMSim(:,tnd2d2post_index);
                    Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
            end
            if schememesh(n) <= d1 + d2
                Yd1d2 = phi(schememesh(n)-d1-d2);
            elseif schememesh(n) > d1 + d2
                if isempty(tnd1d2pre_index); Yd1d2pre = phi(0); else; Yd1d2pre = MMSim(:,tnd1d2pre_index); end; Yd1d2post = MMSim(:,tnd1d2post_index);
                    Yd1d2 = Yd1d2pre + (Yd1d2post-Yd1d2pre)*((tn-d1-d2)-tnd1d2pre)/hsize;
            end
            Yd2d1 = Yd1d2;
            %a,b1,b2,...
            a = B0(Yn,Yd1,Yd2);
            b1 = B1(Yn,Yd1,Yd2); b1x = B1x(Yn,Yd1,Yd2);
            b2 = B2(Yn,Yd1,Yd2); b2x = B2x(Yn,Yd1,Yd2);
            I10Sim = DelW1(1,n)*hsize/2; I01Sim = DelW1(1,n)*hsize - I10Sim; I20Sim = DelW2(1,n)*hsize/2; I02Sim = DelW2(1,n)*hsize - I20Sim;
            expOmega2 = expm( A0hat*hsize + A1*DelW1(1,n) + A2*DelW2(1,n) + 1/2*( LieA0A1*(I10Sim-I01Sim)+LieA0A2*(I20Sim-I02Sim)+LieA1A2*(I21Sim-I12Sim) ) );
            Hx1b1 = (b1x)*(A1*Yn+b1)-A1*b1; % This 0 is actually the Jacobian, d(b1)/dx.   j = 1, l = 1
            Hx1b2 = (b1x)*(A2*Yn+b2)-A2*b1; %  j = 1, l = 2
            Hx2b1 = (b2x)*(A1*Yn+b1)-A1*b2; %  j = 2, l = 1
            Hx2b2 = (b2x)*(A2*Yn+b2)-A2*b2; %  j = 2, l = 1
            b1d1 = B1(Yd1,Yd1d1,Yd2d1); b2d1 = B2(Yd1,Yd1d1,Yd2d1);
            Hxd_one1b1 = B1xd1(Yn,Yd1,Yd2)*(A1*Yd1+b1d1);  % j = 1; l = 1;
            Hxd_one1b2 = B1xd1(Yn,Yd1,Yd2)*(A2*Yd1+b2d1);  % j = 1; l = 2;
            Hxd_one2b1 = B2xd1(Yn,Yd1,Yd2)*(A1*Yd1+b1d1);  % j = 2; l = 1;
            Hxd_one2b2 = B2xd1(Yn,Yd1,Yd2)*(A2*Yd1+b2d1);  % j = 2; l = 2;
            b1d2 = B1(Yd2,Yd1d2,Yd2d2); b2d2 = B2(Yd2,Yd1d2,Yd2d2);
            Hxd_two1b1 = B1xd2(Yn,Yd1,Yd2)*(A1*Yd2+b1d2);  % j = 1; l = 1;
            Hxd_two1b2 = B1xd2(Yn,Yd1,Yd2)*(A2*Yd2+b2d2);  % j = 1; l = 2;
            Hxd_two2b1 = B2xd2(Yn,Yd1,Yd2)*(A1*Yd2+b1d2);  % j = 2; l = 1;
            Hxd_two2b2 = B2xd2(Yn,Yd1,Yd2)*(A2*Yd2+b2d2);  % j = 2; l = 2;
            MMSim(:,n+1) = expOmega2*(Yn + (a - A1*b1 - A2*b2)*hsize + b1*DelW1(1,n) + b2*DelW2(1,n) + Hx1b1*I11Sim + Hx1b2*I21Sim + Hx2b1*I12Sim + Hx2b2*I22Sim + Hxd_one1b1*I11d1Sim + Hxd_one2b1*I12d1Sim + Hxd_one1b2*I21d1Sim + Hxd_one2b2*I22d1Sim  +  Hxd_two1b1*I11d2Sim + Hxd_two2b1*I12d2Sim + Hxd_two1b2*I21d2Sim + Hxd_two2b2*I22d2Sim );

            % Refined Schemes:

                % MilRef:      
                Yn = MilRef(:,n);
                if schememesh(n) <= d1 && schememesh(n) <= d2
                    Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && schememesh(n) <= d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilRef(:,tnd1pre_index); end; Yd1post = MilRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) > 2*d1 && schememesh(n) <= d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilRef(:,tnd1pre_index); end; Yd1post = MilRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MilRef(:,tnd1d1pre_index); end; Yd1d1post = MilRef(:,tnd1d1post_index);
                        Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                    Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) <= d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilRef(:,tnd2pre_index); end; Yd2post = MilRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilRef(:,tnd1pre_index); end; Yd1post = MilRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilRef(:,tnd2pre_index); end; Yd2post = MilRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) > 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilRef(:,tnd1pre_index); end; Yd1post = MilRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MilRef(:,tnd1d1pre_index); end; Yd1d1post = MilRef(:,tnd1d1post_index);
                        Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilRef(:,tnd2pre_index); end; Yd2post = MilRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) <= d1 && 2*d2 < schememesh(n)
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilRef(:,tnd2pre_index); end; Yd2post = MilRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MilRef(:,tnd2d2pre_index); end; Yd2d2post = MilRef(:,tnd2d2post_index);
                        Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                    Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1);
                elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && 2*d2 < schememesh(n)
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilRef(:,tnd1pre_index); end; Yd1post = MilRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilRef(:,tnd2pre_index); end; Yd2post = MilRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MilRef(:,tnd2d2pre_index); end; Yd2d2post = MilRef(:,tnd2d2post_index);
                        Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                    Yd1d1 = phi(schememesh(n)-d1-d1);
                elseif schememesh(n) > 2*d1 && 2*d2 < schememesh(n)
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MilRef(:,tnd1pre_index); end; Yd1post = MilRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MilRef(:,tnd1d1pre_index); end; Yd1d1post = MilRef(:,tnd1d1post_index);
                        Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MilRef(:,tnd2pre_index); end; Yd2post = MilRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MilRef(:,tnd2d2pre_index); end; Yd2d2post = MilRef(:,tnd2d2post_index);
                        Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                end
                if schememesh(n) <= d1 + d2
                    Yd1d2 = phi(schememesh(n)-d1-d2);
                elseif schememesh(n) > d1 + d2
                    if isempty(tnd1d2pre_index); Yd1d2pre = phi(0); else; Yd1d2pre = MilRef(:,tnd1d2pre_index); end; Yd1d2post = MilRef(:,tnd1d2post_index);
                        Yd1d2 = Yd1d2pre + (Yd1d2post-Yd1d2pre)*((tn-d1-d2)-tnd1d2pre)/hsize;
                end
                Yd2d1 = Yd1d2;
                %a,b1,b2,...
                a = A0*Yn + B0(Yn,Yd1,Yd2);
                b1 = A1*Yn + B1(Yn,Yd1,Yd2); b1x = A1+B1x(Yn,Yd1,Yd2);
                b2 = A2*Yn + B2(Yn,Yd1,Yd2); b2x = A2+B2x(Yn,Yd1,Yd2);
                I11Ref = (DelW1(1,n)^2-hsize)/2; I22Ref = (DelW2(1,n)^2-hsize)/2;
                    % I12ref calculation:
                    I12Ref = dW1(1,i1f)*dW2(1,i1f)/2;
                    for j=i1f+1:i2f-1
                        I12Ref = I12Ref + dW1(1,j)*dW2(1,j)/2 + dW2(1,j)*(W1(1,j)-W1(1,i1f));
                    end
                    I21Ref = DelW1(1,n)*DelW2(1,n) - I12Ref;
                Lx1b1 = (b1x)*b1;   Lx2b2 = (b2x)*b2;   Lx1b2 = (b2x)*b1;   Lx2b1 = (b1x)*b2;
                b1xd1 = B1xd1(Yn,Yd1,Yd2); b2xd1 = B2xd1(Yn,Yd1,Yd2); b1d1 = A1*Yd1 + B1(Yd1,Yd1d1,Yd2d1); b2d1 = A2*Yd1 + B2(Yd1,Yd1d1,Yd2d1);
                Lxd_one1b1 = (b1xd1)*b1d1; Lxd_one2b2 = (b2xd1)*b2d1;   Lxd_one1b2 = (b2xd1)*b1d1;   Lxd_one2b1 = (b1xd1)*b2d1;
                % Iijd1Ref calculations:
                if schememesh(n) < d1
                    I11d1Ref = 0; I12d1Ref = 0; I21d1Ref = 0; I22d1Ref = 0;
                else
                    % I11d1Ref:
                        I11d1Ref = dW1d1(1,i1f)*dW1(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I11d1Ref = I11d1Ref + dW1d1(1,j)*dW1(1,j)/2 + dW1(1,j)*(W1d1(1,j)-W1d1(1,i1f));
                        end
                    % I12d1Ref:
                        I12d1Ref = dW1d1(1,i1f)*dW2(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I12d1Ref = I12d1Ref + dW1d1(1,j)*dW2(1,j)/2 + dW2(1,j)*(W1d1(1,j)-W1d1(1,i1f));
                        end
                    % I21d1Ref:
                        I21d1Ref = dW1(1,i1f)*dW2d1(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I21d1Ref = I21d1Ref + dW1(1,j)*dW2d1(1,j)/2 + dW2d1(1,j)*(W1(1,j)-W1(1,i1f));
                        end
                    % I22d1Ref:
                        I22d1Ref = dW2(1,i1f)*dW2d1(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I22d1Ref = I22d1Ref + dW2(1,j)*dW2d1(1,j)/2 + dW2d1(1,j)*(W2(1,j)-W2(1,i1f));
                        end
                    I21d1Ref = DelW1(1,n)*DelW2d1(1,n)-I21d1Ref;
                    I22d1Ref = DelW2(1,n)*DelW2d1(1,n)-I22d1Ref;
                end
                b1xd2 = B1xd2(Yn,Yd1,Yd2); b2xd2 = B2xd2(Yn,Yd1,Yd2); b1d2 = A1*Yd2 + B1(Yd2,Yd1d2,Yd2d2); b2d2 = A2*Yd2 + B2(Yd2,Yd1d2,Yd2d2);
                Lxd_two1b1 = (b1xd2)*b1d2; Lxd_two2b2 = (b2xd2)*b2d2;   Lxd_two1b2 = (b2xd2)*b1d2;   Lxd_two2b1 = (b1xd2)*b2d2;
                % Iijd2Ref calculations:
                if schememesh(n) < d2
                    I11d2Ref = 0; I12d2Ref = 0; I21d2Ref = 0; I22d2Ref = 0;
                else
                    % I11d2Ref:
                        I11d2Ref = dW1d2(1,i1f)*dW1(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I11d2Ref = I11d2Ref + dW1d2(1,j)*dW1(1,j)/2 + dW1(1,j)*(W1d2(1,j)-W1d2(1,i1f));
                        end
                    % I12d2Ref:
                        I12d2Ref = dW1d2(1,i1f)*dW2(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I12d2Ref = I12d2Ref + dW1d2(1,j)*dW2(1,j)/2 + dW2(1,j)*(W1d2(1,j)-W1d2(1,i1f));
                        end
                    % I21d2Ref:
                        I21d2Ref = dW1(1,i1f)*dW2d2(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I21d2Ref = I21d2Ref + dW1(1,j)*dW2d2(1,j)/2 + dW2d2(1,j)*(W1(1,j)-W1(1,i1f));
                        end
                    % I22d2Ref:
                        I22d2Ref = dW2(1,i1f)*dW2d2(1,i1f)/2;
                        for j=i1f+1:i2f-1
                            I22d2Ref = I22d2Ref + dW2(1,j)*dW2d2(1,j)/2 + dW2d2(1,j)*(W2(1,j)-W2(1,i1f));
                        end
                    I21d2Ref = DelW1(1,n)*DelW2d2(1,n)-I21d2Ref;
                    I22d2Ref = DelW2(1,n)*DelW2d2(1,n)-I22d2Ref;
                end
                MilRef(:,n+1) = Yn + a*hsize + b1*DelW1(1,n) + b2*DelW2(1,n) + Lx1b1*I11Ref + Lx2b1*I21Ref + Lx1b2*I12Ref + Lx2b2*I22Ref + Lxd_one1b1*I11d1Ref + Lxd_one1b2*I12d1Ref + Lxd_one2b1*I21d1Ref + Lxd_one2b2*I22d1Ref  +  Lxd_two1b1*I11d2Ref + Lxd_two1b2*I12d2Ref + Lxd_two2b1*I21d2Ref + Lxd_two2b2*I22d2Ref    ;

                % MMRef:
                Yn = MMRef(:,n);
                if schememesh(n) <= d1 && schememesh(n) <= d2
                    Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && schememesh(n) <= d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMRef(:,tnd1pre_index); end; Yd1post = MMRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    Yd1d1 = phi(schememesh(n)-d1-d1); Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) > 2*d1 && schememesh(n) <= d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMRef(:,tnd1pre_index); end; Yd1post = MMRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MMRef(:,tnd1d1pre_index); end; Yd1d1post = MMRef(:,tnd1d1post_index);
                        Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                    Yd2 = phi(schememesh(n)-d2); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) <= d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMRef(:,tnd2pre_index); end; Yd2post = MMRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMRef(:,tnd1pre_index); end; Yd1post = MMRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMRef(:,tnd2pre_index); end; Yd2post = MMRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    Yd1d1 = phi(schememesh(n)-d1-d1); Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) > 2*d1 && d2 < schememesh(n) && schememesh(n) <= 2*d2
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMRef(:,tnd1pre_index); end; Yd1post = MMRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MMRef(:,tnd1d1pre_index); end; Yd1d1post = MMRef(:,tnd1d1post_index);
                        Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMRef(:,tnd2pre_index); end; Yd2post = MMRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    Yd2d2 = phi(schememesh(n)-d2-d2);
                elseif schememesh(n) <= d1 && 2*d2 < schememesh(n)
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMRef(:,tnd2pre_index); end; Yd2post = MMRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MMRef(:,tnd2d2pre_index); end; Yd2d2post = MMRef(:,tnd2d2post_index);
                        Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                    Yd1 = phi(schememesh(n)-d1); Yd1d1 = phi(schememesh(n)-d1-d1);
                elseif d1 < schememesh(n) && schememesh(n) <= 2*d1 && 2*d2 < schememesh(n)
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMRef(:,tnd1pre_index); end; Yd1post = MMRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMRef(:,tnd2pre_index); end; Yd2post = MMRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MMRef(:,tnd2d2pre_index); end; Yd2d2post = MMRef(:,tnd2d2post_index);
                        Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                    Yd1d1 = phi(schememesh(n)-d1-d1);
                elseif schememesh(n) > 2*d1 && 2*d2 < schememesh(n)
                    if isempty(tnd1pre_index); Yd1pre = phi(0); else; Yd1pre = MMRef(:,tnd1pre_index); end; Yd1post = MMRef(:,tnd1post_index);
                        Yd1 = Yd1pre + (Yd1post-Yd1pre)*((tn-d1)-tnd1pre)/hsize;
                    if isempty(tnd1d1pre_index); Yd1d1pre = phi(0); else; Yd1d1pre = MMRef(:,tnd1d1pre_index); end; Yd1d1post = MMRef(:,tnd1d1post_index);
                        Yd1d1 = Yd1d1pre + (Yd1d1post-Yd1d1pre)*((tn-d1-d1)-tnd1d1pre)/hsize;
                    if isempty(tnd2pre_index); Yd2pre = phi(0); else; Yd2pre = MMRef(:,tnd2pre_index); end; Yd2post = MMRef(:,tnd2post_index);
                        Yd2 = Yd2pre + (Yd2post-Yd2pre)*((tn-d2)-tnd2pre)/hsize;
                    if isempty(tnd2d2pre_index); Yd2d2pre = phi(0); else; Yd2d2pre = MMRef(:,tnd2d2pre_index); end; Yd2d2post = MMRef(:,tnd2d2post_index);
                        Yd2d2 = Yd2d2pre + (Yd2d2post-Yd2d2pre)*((tn-d2-d2)-tnd2d2pre)/hsize;
                end
                if schememesh(n) <= d1 + d2
                    Yd1d2 = phi(schememesh(n)-d1-d2);
                elseif schememesh(n) > d1 + d2
                    if isempty(tnd1d2pre_index); Yd1d2pre = phi(0); else; Yd1d2pre = MMRef(:,tnd1d2pre_index); end; Yd1d2post = MMRef(:,tnd1d2post_index);
                        Yd1d2 = Yd1d2pre + (Yd1d2post-Yd1d2pre)*((tn-d1-d2)-tnd1d2pre)/hsize;
                end
                Yd2d1 = Yd1d2;
                %a,b1,b2,...
                a = B0(Yn,Yd1,Yd2);
                b1 = B1(Yn,Yd1,Yd2); b1x = B1x(Yn,Yd1,Yd2);
                b2 = B2(Yn,Yd1,Yd2); b2x = B2x(Yn,Yd1,Yd2);
                I10Ref = dW1(1,i1f)*(t(i1f+1)-t(i1f))/2;
                    for j=i1f+1:i2f-1
                        I10Ref = I10Ref + dW1(1,j)*(t(i1f+1)-t(i1f))/2 + refdt*(W1(1,j)-W1(1,i1f));
                    end
                I01Ref = DelW1(1,n)*hsize - I10Ref;
                I20Ref = dW2(1,i1f)*(t(i1f+1)-t(i1f))/2;
                    for j=i1f+1:i2f-1
                        I20Ref = I20Ref + dW2(1,j)*(t(i1f+1)-t(i1f))/2 + refdt*(W2(1,j)-W2(1,i1f));
                    end
                I02Ref = DelW2(1,n)*hsize - I20Ref;
                expOmega2 = expm( A0hat*hsize + A1*DelW1(1,n) + A2*DelW2(1,n) + 1/2*( LieA0A1*(I10Ref-I01Ref)+LieA0A2*(I20Ref-I02Ref)+LieA1A2*(I21Ref-I12Ref) ) );
                Hx1b1 = (b1x)*(A1*Yn+b1)-A1*b1; % This 0 is actually the Jacobian, d(b1)/dx.   j = 1, l = 1
                Hx1b2 = (b1x)*(A2*Yn+b2)-A2*b1; %  j = 1, l = 2
                Hx2b1 = (b2x)*(A1*Yn+b1)-A1*b2; %  j = 2, l = 1
                Hx2b2 = (b2x)*(A2*Yn+b2)-A2*b2; %  j = 2, l = 1
                b1d1 = B1(Yd1,Yd1d1,Yd2d1); b2d1 = B2(Yd1,Yd1d1,Yd2d1);
                Hxd_one1b1 = B1xd1(Yn,Yd1,Yd2)*(A1*Yd1+b1d1);  % j = 1; l = 1;
                Hxd_one1b2 = B1xd1(Yn,Yd1,Yd2)*(A2*Yd1+b2d1);  % j = 1; l = 2;
                Hxd_one2b1 = B2xd1(Yn,Yd1,Yd2)*(A1*Yd1+b1d1);  % j = 2; l = 1;
                Hxd_one2b2 = B2xd1(Yn,Yd1,Yd2)*(A2*Yd1+b2d1);  % j = 2; l = 2;

                b1d2 = B1(Yd2,Yd1d2,Yd2d2); b2d2 = B2(Yd2,Yd1d2,Yd2d2);
                Hxd_two1b1 = B1xd2(Yn,Yd1,Yd2)*(A1*Yd2+b1d2);  % j = 1; l = 1;
                Hxd_two1b2 = B1xd2(Yn,Yd1,Yd2)*(A2*Yd2+b2d2);  % j = 1; l = 2;
                Hxd_two2b1 = B2xd2(Yn,Yd1,Yd2)*(A1*Yd2+b1d2);  % j = 2; l = 1;
                Hxd_two2b2 = B2xd2(Yn,Yd1,Yd2)*(A2*Yd2+b2d2);  % j = 2; l = 2;

                MMRef(:,n+1) = expOmega2*(Yn + (a - A1*b1 - A2*b2)*hsize + b1*DelW1(1,n) + b2*DelW2(1,n) + Hx1b1*I11Ref + Hx2b1*I12Ref + Hx1b2*I21Ref + Hx2b2*I22Ref + Hxd_one1b1*I11d1Ref + Hxd_one2b1*I12d1Ref + Hxd_one1b2*I21d1Ref + Hxd_one2b2*I22d1Ref  +  Hxd_two1b1*I11d2Ref + Hxd_two2b1*I12d2Ref + Hxd_two1b2*I21d2Ref + Hxd_two2b2*I22d2Ref  );
                
            
        end

%% SAVING VALUES:        

        if trial == samplepath && hvalues(h) == samplehvalue
            sampletime = schememesh;
            sampleEM = EM(samplecomponent,:);
            sampleMilSim = MilSim(samplecomponent,:);
            sampleMEM = MEM(samplecomponent,:);
            sampleMMSim = MMSim(samplecomponent,:);
            sampleMilRef = MilRef(samplecomponent,:);
            sampleMMRef = MMRef(samplecomponent,:);
        end

        for ti = 1:NO
%             if isempty(find(schememesh==schememesh(schememeshindexobservations(ti)),1))
%                 ipre = find(schememesh<=schememesh(schememeshindexobservations(ti)),1,'last');
%                 ipost = find(schememesh>schememesh(tindexobservations(ti)),1);
%                 EMvalues(:,ti,trial,h) = EM(:,ipre) + (EM(:,ipost)-EM(:,ipre))*(t(tindexobservations(ti))-schememesh(ipre))/Delt;
%                 MilSimvalues(:,ti,trial,h) = MilSim(:,ipre) + (MilSim(:,ipost)-MilSim(:,ipre))*(t(tindexobservations(ti))-schememesh(ipre))/Delt;
%                 MEMvalues(:,ti,trial,h) = MEM(:,ipre) + (MEM(:,ipost)-MEM(:,ipre))*(t(tindexobservations(ti))-schememesh(ipre))/Delt;
%                 MMSimvalues(:,ti,trial,h) = MMSim(:,ipre) + (MMSim(:,ipost)-MMSim(:,ipre))*(t(tindexobservations(ti))-schememesh(ipre))/Delt;
%                 MilRefvalues(:,ti,trial,h) = MilRef(:,ipre) + (MilRef(:,ipost)-MilRef(:,ipre))*(t(tindexobservations(ti))-schememesh(ipre))/Delt;
%                 MMRefvalues(:,ti,trial,h) = MMRef(:,ipre) + (MMRef(:,ipost)-MMRef(:,ipre))*(t(tindexobservations(ti))-schememesh(ipre))/Delt;
%             else
                i3 = find(schememesh==schememesh(schememeshindexobservations(ti)),1);
                EMvalues(:,ti,trial,h) = EM(:,i3);
                MilSimvalues(:,ti,trial,h) = MilSim(:,i3);
                MEMvalues(:,ti,trial,h) = MEM(:,i3);
                MMSimvalues(:,ti,trial,h) = MMSim(:,i3);
                MilRefvalues(:,ti,trial,h) = MilRef(:,i3);
                MMRefvalues(:,ti,trial,h) = MMRef(:,i3);
%             end
        end

        for ti = 1:NO
            ErrorEM(ti,trial,h) = norm(EMvalues(:,ti,trial,h)-Xrefvalues(:,ti,trial));
            ErrorMilSim(ti,trial,h) = norm(MilSimvalues(:,ti,trial,h)-Xrefvalues(:,ti,trial));
            ErrorMEM(ti,trial,h) = norm(MEMvalues(:,ti,trial,h)-Xrefvalues(:,ti,trial));
            ErrorMMSim(ti,trial,h) = norm(MMSimvalues(:,ti,trial,h)-Xrefvalues(:,ti,trial));
            ErrorMilRef(ti,trial,h) = norm(MilRefvalues(:,ti,trial,h)-Xrefvalues(:,ti,trial));
            ErrorMMRef(ti,trial,h) = norm(MMRefvalues(:,ti,trial,h)-Xrefvalues(:,ti,trial));
        end

    end

end

%% ERROR CALCULATIONS

for ti = 1:NO
    for h = 1:length(hvalues)
        MSErrorEM(ti,h) = sqrt( mean( squeeze( ErrorEM(ti,:,h).^2 ) ) );
        MSErrorMilSim(ti,h) = sqrt( mean( squeeze( ErrorMilSim(ti,:,h).^2 ) ) );
        MSErrorMEM(ti,h) = sqrt( mean( squeeze( ErrorMEM(ti,:,h).^2 ) ) );
        MSErrorMMSim(ti,h) = sqrt( mean( squeeze( ErrorMMSim(ti,:,h).^2 ) ) );
        MSErrorMilRef(ti,h) = sqrt( mean( squeeze( ErrorMilRef(ti,:,h).^2 ) ) );
        MSErrorMMRef(ti,h) = sqrt( mean( squeeze( ErrorMMRef(ti,:,h).^2 ) ) );
    end
end

%% TRAJECTORY PLOT:

figure(6+simulation)
hold on
pretime = -1:refdt:0; history = zeros(1,length(pretime)); for j = 1:length(pretime); phi2 = phi(pretime(j)); history(j) = phi2(samplecomponent); end
plot([pretime,t],[history,sampleXref],'k')
plot([0,sampletime],[X0(samplecomponent),sampleEM],'r-o')
plot([0,sampletime],[X0(samplecomponent),sampleMilSim],'-o','color',[1/3,1/3,1])
     plot([0,sampletime],[X0(samplecomponent),sampleMilRef],'-square','color',[1/3,1/3,1])
plot([0,sampletime],[X0(samplecomponent),sampleMEM],'--*','color',[1,1/2,0])
plot([0,sampletime],[X0(samplecomponent),sampleMMSim],'--*','color',[0,2/3,0])
     plot([0,sampletime],[X0(samplecomponent),sampleMMRef],'--square','MarkerSize',10,'color',[0,2/3,0])
     legend('Reference','EM','Milstein (Simple)','Milstein (Refined)','MEM','MM (Simple)','MM (Refined)','Interpreter','latex','FontSize',12,'location','Northwest')
xlabel('$t$','Interpreter','latex','FontSize',12)
ylabel('$X_1(t)$','Interpreter','latex','FontSize',12)
title('Sample Trajectory')
axis([-1,T,min(sampleXref)-0.15,max(sampleXref)+0.15])
set(gcf,'position',[200,75,1200,500])


%% ERROR GRAPHS:

if showerrorplots == 1
% Error graph (t=Ti):
for Ti = 1:NO
    figure(simulation)
    hold on
    %referenceorderhalf = max( [MSErrorEM(Ti,1),MSErrorMEM(Ti,1),MSErrorMilSim(Ti,1)] )*sqrt(hvalues); referenceorderone = min( [MSErrorEM(Ti,1),MSErrorMEM(Ti,1),MSErrorMilSim(Ti,1),MSErrorMMSim(Ti,1)] )*hvalues;
    referenceorderhalf = sqrt(hvalues); referenceorderone = hvalues;
    if simulation == 1      
        plot(log(hvalues)/log(2),log10(referenceorderhalf)-0.5,'m-o','Linewidth',1)
        plot(log(hvalues)/log(2),log10(referenceorderone)-1.25,'m-square','Linewidth',1)
    elseif simulation == 2
        plot(log(hvalues)/log(2),log10(referenceorderhalf)-0.5,'m-o','Linewidth',1)
        plot(log(hvalues)/log(2),log10(referenceorderone)-1,'m-square','Linewidth',1)
    elseif simulation == 3
        plot(log(hvalues)/log(2),log10(referenceorderhalf)-0.25,'m-o','Linewidth',1)
        plot(log(hvalues)/log(2),log10(referenceorderone)-1,'m-square','Linewidth',1)
    end
    plot(log(hvalues)/log(2),log(MSErrorEM(Ti,:))/log(10),'r-o','Linewidth',1)
    plot(log(hvalues)/log(2),log(MSErrorMilSim(Ti,:))/log(10),'-o','color',[1/3,1/3,1],'Linewidth',1)
        plot(log(hvalues)/log(2),log(MSErrorMilRef(Ti,:))/log(10),'-square','color',[1/3,1/3,1],'Linewidth',1)
    plot(log(hvalues)/log(2),log(MSErrorMEM(Ti,:))/log(10),'--*','color',[1,1/2,0],'Linewidth',1)
    plot(log(hvalues)/log(2),log(MSErrorMMSim(Ti,:))/log(10),'--*','color',[0,2/3,0],'Linewidth',1)
        plot(log(hvalues)/log(2),log(MSErrorMMRef(Ti,:))/log(10),'--square','MarkerSize',10,'color',[0,2/3,0],'Linewidth',1)
    xticks(log(dt)/log(2):0)
    xlabel('$\log_2h$','Interpreter','latex','FontSize',16)
    % if mod(T/d1,1) == 0 && mod(T/d2,1) == 0
        legend('Reference $1/2$','Reference $1$','EM LI','Milstein LI (Simple)','Milstein LI (Refined)','MEM LI','MM LI (Simple)','MM LI (Refined)','Interpreter','latex','FontSize',14,'location','Southeast','NumColumns',1)
    % else
    %     legend('EM','Milstein (Simple)','Milstein (Refined)','MEM','MM (Simple)','MM (Refined)','Interpreter','latex','FontSize',14,'location','Southeast','NumColumns',1)
    % end
        %axis([log(dt)/log(2),log(1)/log(2),axisbounds(1,Ti)-0.5,axisbounds(2,Ti)+0.25])
    xlim([log(dt)/log(2)-0.1,max(log(hvalues)/log(2))+0.1])
    ylabel('$\log_{10}\mathrm{MSE}(4)$','Interpreter','latex','FontSize',16)
    if t(tindexobservations(Ti)) == round(t(tindexobservations(Ti)),4)
        %title(['MS Strong Errors ','($t=$ ',num2str(t(tindexobservations(Ti))),')'],'Interpreter','latex','FontSize',16)
    elseif t(tindexobservations(Ti)) ~= round(t(tindexobservations(Ti)),4)
        %title(['MS Strong Errors ','($t\approx$ ',num2str(t(tindexobservations(Ti))),')'],'Interpreter','latex','FontSize',16)
    end
    %axis([log(dt)/log(2)-0.1,max(log(hvalues)/log(2))+0.1,-5.52,-0.95])
    %axis([log(dt)/log(2)-0.1,max(log(hvalues)/log(2))+0.1,-5,-0.5])
    set(gcf,'position',[200,75,700,500])
    grid on
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)
    b = get(gca,'YTickLabel');
    set(gca,'YTickLabel',b,'fontsize',12)
    % EXTRA SAVING:
    % COMMENT THIS OUT AFTER USE:
    if saveplots == 1 && Ti == NO && runs == runstosave && simulation == 1
        saveas(gcf,'ZExample1LIA.eps','epsc')
    elseif saveplots == 1 && Ti == NO && runs == runstosave && simulation == 2
        legend off
        saveas(gcf,'ZExample1LIB.eps','epsc')
    elseif saveplots == 1 && Ti == NO && runs == runstosave && simulation == 3
        legend off
        saveas(gcf,'ZExample1LIC.eps','epsc')
    elseif saveplots == 1 && Ti == NO && runs == runstosave && simulation == 4
        legend off
        saveas(gcf,'ZExample1LID.eps','epsc')
    elseif saveplots == 1 && Ti == NO && runs == runstosave && simulation == 5
        legend off
        saveas(gcf,'ZExample1LIE.eps','epsc')
    elseif saveplots == 1 && Ti == NO && runs == runstosave && simulation == 6
        legend off
        saveas(gcf,'ZExample1LIF.eps','epsc')
    end
end
end

end

disp(['Total simulation time is ',num2str(round(toc(simulationtime)/60,4)),' minutes.'])