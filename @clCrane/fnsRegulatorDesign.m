function [t, x, u, idxExit] = fnsRegulatorDesign(A, B, t, x0, sIdx, v0, v1, v2, ...
    oG, oq, qG, q0, poles)

% Useless comment for useless change

% This experimental code simulates the solution using regulator-simplex
% design. Outside code will have to manage initial conditions and piecing
% together the solution from consecutive simplices

% Inputs:
% A, B: system matrices
% t:    time over which to solve trajectory for
% x0:   initial state for trajectory (4x1 vector)
% sIdx: the index of the current simplex (for plotting annotation)
% v0, v1, v2: (x,y) coordinates of the simplex, F0 is the target facet
% oG:   Indices of unknowns s1-s6 that will be guessed
% oq:   Indices of unknowns s1-s6 that will be solved for
% qG:   Corresponding guessed values for oG
% q0:   Initial values for fsolve for oq
% poles:    usual pole placement component, for tracking convergence

C = [1 0 0 0; 0 0 1 0]; % system output
Cr = eye(2); % exosystem output

dispDump = 0;

fprintf('Simplex %d\n', sIdx);

%% Design an exosystem with the desired trajectories on the simplex

h1 = findOutwardNormalVector(v0, v2, v1);
h2 = findOutwardNormalVector(v0, v1, v2);
xi = obtainFlowVector(v0, v1, v2);

if dispDump
    % Write in terms of symbolic variables
    syms s1 s2 s3 s4 s5 s6
    Ss = [s1 s2; s3 s4];
    ss = [s5; s6];
    y0s = Ss * v0 + ss;
    y1s = Ss * v1 + ss;
    y2s = Ss * v2 + ss;
    % Invariance conditions
    evalIneq(h1, y0s, 'h1 . y0', '<=', 1);
    evalIneq(h2, y0s, 'h2 . y0', '<=', 1);
    evalIneq(h2, y1s, 'h2 . y1', '<=', 1);
    evalIneq(h1, y2s, 'h1 . y2', '<=', 1);
    % Flow conditions
    evalIneq(xi, y0s, 'xi . y0', '<', 1);
    evalIneq(xi, y1s, 'xi . y1', '<', 1);
    evalIneq(xi, y2s, 'xi . y2', '<', 1);
end

% Compute exosystem that solves FBI equations by construction
[S, s] = computeExosystembyFBI(A, B, x0, oG, oq, qG, q0)
%[S, s] = computeExosystembyFBI2(A, B, x0, q0)

% Initial condition for exosystem
w0 = [x0(1); x0(3)];

% Compute the controller that satisfies solves the (ATP)
[PI, GAMMA, p, gm] = findController(S, s, A, B);

% Sanity checks
vS = [S(1,1)*w0(1) + S(1,2)*w0(2) + s(1); ...
    S(2,1)*w0(1) + S(2,2)*w0(2) + s(2)];
fprintf('Velocity difference norm = %f\n', norm(vS - x0([2,4])));
fprintf('FBI Equations max norm from 0: %f\n', ...
    FBISolvedNorm(A, B, C, S, s, Cr, PI, GAMMA, p, gm));

y0 = S * v0 + s;
y1 = S * v1 + s;
y2 = S * v2 + s;

% Invariance conditions
evalIneq(h1, y0, 'h1 . y0', '<=', 1);
evalIneq(h2, y0, 'h2 . y0', '<=', 1);
evalIneq(h2, y1, 'h2 . y1', '<=', 1);
evalIneq(h1, y2, 'h1 . y2', '<=', 1);

% Flow conditions
% Technically not treated correctly, because can any xi, not just a fixed
% one; equivalently, checking for equilibrium points inside the simplex
% will tell whether flow conditions hold
evalIneq(xi, y0, 'xi . y0', '<', 1);
evalIneq(xi, y1, 'xi . y1', '<', 1);
evalIneq(xi, y2, 'xi . y2', '<', 1);

%% Simulate the trajectory for the (A,B) system

K = -place(A, B, poles);

% Simulate the composite system with the exosystem
Asys = [A + B*K, B*(GAMMA-K*PI); zeros(2,4), S];
asys = [B*(gm - K*p); s];
[~, x] = ode45(@(t, x) regulatorODE(t, x, Asys, asys), t, ...
    [x0; w0], []);

%% Determine the state vector when it exited the simplex

idxExit = [];
[xint, yint] = polyxpoly([v1(1); v2(1)], [v1(2); v2(2)], x(:,1), x(:,3));
if ~isempty(xint)   
    exitPts = zeros(size(xint, 1), 1);
    for j = 1:size(xint, 1)
        distBest = 1e6;
        for i = 1:size(x, 1)
            dist = norm((x(i, [1,3]))' - [xint(j); yint(j)]);
            if dist < distBest
                distBest = dist;
                exitPts(j) = i;
            end
        end       
    end
    idxExit = min(exitPts);
end

%% Create plots

scrSz = get(0,'ScreenSize');
hDraw = figure('units', 'normalized', 'outerposition', [0 0 scrSz(4)/scrSz(3) 1]);

% Domain over which to solve the vector field
maxx = max([v0(1), v1(1), v2(1)]);
minx = min([v0(1), v1(1), v2(1)]);
deltax = 0.3 * (maxx - minx);
maxx = maxx + deltax;
minx = minx - deltax;
maxy = max([v0(2), v1(2), v2(2)]);
miny = min([v0(2), v1(2), v2(2)]);
deltay = 0.3 * (maxy - miny);
maxy = maxy + deltay;
miny = miny - deltay;

% Determine equilibrium point(s)
switch(rank(S))
    case 0
        % S = 0, nothing moves
    case 1
        % Eqb line
        % Several expressions are equivalent
        wEqb(1,:) = minx:0.1:maxx;
        wEqb(2,:) = -1/S(2,2) * (S(2,1) * wEqb(1,:) + s(2));
    case 2      
        wEqb = S \ (-s); % equilibrium point     
        
        % Expand domain
        maxx = max(maxx, wEqb(1) + deltax);
        maxy = max(maxy, wEqb(2) + deltay);
        minx = min(minx, wEqb(1) - deltax);
        miny = min(miny, wEqb(2) - deltay);
end

% Plot the exosystem vector field over the simplex
subplot(2, 2, 1);
hold all;
[X,Y] = meshgrid(minx:0.1:maxx, miny:0.1:maxy);
streamslice(X, Y, S(1, 1) * X + S(1, 2) * Y + s(1), S(2, 1) * X + S(2, 2) * Y + s(2));
line([v0(1), v1(1)], [v0(2), v1(2)], 'Color', 'k', 'LineWidth', 2);
line([v2(1), v1(1)], [v2(2), v1(2)], 'Color', [1.0000 0.6471 0], ...
    'LineWidth', 2); % F0
line([v0(1), v2(1)], [v0(2), v2(2)], 'Color', 'k', 'LineWidth', 2);
% quiver(mean([v0(1),v1(1)]), mean([v0(2),v1(2)]), h2(1), h2(2), ...
%     'Color', 'k', 'LineWidth', 1);
% quiver(mean([v0(1),v2(1)]), mean([v0(2),v2(2)]), h1(1), h1(2), ...
%     'Color', 'k', 'LineWidth', 1);
% quiver(v0(1), v0(2), xi(1), xi(2), 'Color', 'c', 'LineWidth', 2);
text(v0(1), v0(2), 'v_0', 'BackgroundColor', [0.4980    1.0000    0.8314]);
text(v1(1), v1(2), 'v_1', 'BackgroundColor', [0.4980    1.0000    0.8314]);
text(v2(1), v2(2), 'v_2', 'BackgroundColor', [0.4980    1.0000    0.8314]);
% text(v0(1)+xi(1), v0(2)+xi(2), '$\xi$', 'BackgroundColor', [.7 .9 .7], ...
%     'Interpreter', 'Latex');

% Plot eqb point(s)
switch(rank(S))
    case 1
        plot(wEqb(1,:), wEqb(2,:), 'm-');
    case 2        
        plot(wEqb(1, :), wEqb(2, :), 'm*');
end

axis equal
xlabel('x_1 (m)');
ylabel('x_3 (rad)');
title(sprintf('Vector field over simplex %d', sIdx));

% Plot the trajectory of exosystem and system over the simplex
subplot(2,2,2);
hold all;

line([v0(1), v1(1)], [v0(2), v1(2)], 'Color', 'k', 'LineWidth', 2);
line([v2(1), v1(1)], [v2(2), v1(2)], 'Color', [1.0000 0.6471 0], ...
    'LineWidth', 2); % F0
line([v0(1), v2(1)], [v0(2), v2(2)], 'Color', 'k', 'LineWidth', 2);

if isempty(idxExit)
    plot(x(:, 5), x(:, 6), 'r--', 'LineWidth', 2);
    plot(x(:, 1), x(:, 3), 'g-', 'LineWidth', 2);
else
    plot(x(1:idxExit, 5), x(1:idxExit, 6), 'r--', 'LineWidth', 2);
    plot(x(1:idxExit, 1), x(1:idxExit, 3), 'g-', 'LineWidth', 2);
end

axis equal
xlabel('x_1 (m)');
ylabel('x_3 (rad)');
title('State Space Tracking');

% Compare trajectories of position
subplot(2,2,3);
plot(t, x(:, 1), 'g-', t, x(:, 5), 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('x_1 (m)');
title('Cart Position vs. Time');
legend('Cart Position', 'Exosystem signal');
grid on;

% Compare trajectories of angle
subplot(2,2,4);
plot(t, x(:, 3), 'g-', t, x(:, 6), 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('x_3 (rad)');
title('Crane Angle vs. Time');
legend('Crane Angle', 'Exosystem signal');
grid on;

% Save plot
set(hDraw,'Units', 'Inches');
pos = get(hDraw, 'Position');
set(hDraw, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', ...
    'PaperSize', [pos(3), pos(4)]);
saveName = sprintf('RegulatorResult_s%d', sIdx);
print(hDraw, '-dpdf', ['PlotDump/' saveName]);
close(hDraw);

%% Last step

% Obtain control input history
nPts = size(x,1);
u = zeros(nPts, 1);
for i = 1:nPts
    u(i) = GAMMA * x(i, 5:6)' + gm + K * (x(i,1:4)' - PI * x(i, 5:6)' - p);
end

x(:,5:6) = []; % do not return the exosystem state

end

function h = findOutwardNormalVector(x1, x2, x3)

% Find the normal to the line segment x1-x2
% The three points x1, x2, x3 define a triangle, so have the normal vector
% point outside

dx = x1 - x2;
h = [0 -1; 1 0] * dx;
h = h / norm(h);
if dot(h, x3 - x1) > 0 % x3 - x2 works too
    h = -h;
end

if dot(h, dx) > 1e-06
    error('Normal vector not found');
end
end

function xi = obtainFlowVector(v0, v1, v2)

% Decided to have this point away from the centroid of the simplex, based
% on the v0 point (exit facet is F0)
V = [v0 v1 v2];
cen = [mean(V(1, :)); mean(V(2, :))];
xi = v0 - cen;
xi = xi / norm(xi);
end

function F = funcFBIConstraints(q, A, B, x0, oG, oq, qG)

% Unpack the variables
T = [oG', oq'; qG', q'];
[~, o] = sort(T(1, :));
T = T(:, o);
s1 = T(2, 1);
s2 = T(2, 2);
s3 = T(2, 3);
s4 = T(2, 4);
s5 = T(2, 5);
s6 = T(2, 6);

F = [B(4)*(s1^2 + s2*s3 - A(2,4)*s3) - B(2)*(s1*s3 + s3*s4 - A(4,4)*s3);...
    B(4)*(s1*s2 + s2*s4 - A(2,3) - A(2,4)*s3) - B(2)*(s2*s3 + s4^2 - A(4,3) - A(4,4)*s4);...
    s5*(B(4)*s1 - B(2)*s3) + s6*(B(4)*(s2 - A(2,4)) - B(2)*(s4 + A(4,4)))];

if size(q, 1) == 5 % velocity matching constraint
    F = [F;...
        x0(2) - s5 - x0(1)*s1 - x0(3)*s2;
        x0(4) - s6 - x0(1)*s3 - x0(3)*s4];
end
end

function F = funcFBIConstraints2(q, A, B, x0)

% Unpack the variables
s1 = q(1);
s2 = q(2);
s3 = q(3);
s4 = q(4);
s5 = q(5);
s6 = q(6);

F = [B(4)*(s1^2 + s2*s3 - A(2,4)*s3) - B(2)*(s1*s3 + s3*s4 - A(4,4)*s3);...
    B(4)*(s1*s2 + s2*s4 - A(2,3) - A(2,4)*s3) - B(2)*(s2*s3 + s4^2 - A(4,3) - A(4,4)*s4);...
    s5*(B(4)*s1 - B(2)*s3) + s6*(B(4)*(s2 - A(2,4)) - B(2)*(s4 + A(4,4)));...
    x0(2) - s5 - x0(1)*s1 - x0(3)*s2; ...
    x0(4) - s6 - x0(1)*s3 - x0(3)*s4; ...
    s1*s4 - s2*s3 - 1]; % determinant not zero 

end

function [S, s] = computeExosystembyFBI(A, B, x0, oG, oq, qG, q0)

idxOrder = sort([oG; oq]);
if ~isequal(idxOrder', 1:6)
    error('Error setting up FBI Constraint solver');
end

% Solve the system of non-linear equations/constraints
[qSol, ~] = fsolve(@(q) funcFBIConstraints(q, A, B, x0, oG, oq, qG), ...
    q0, optimset('Display', 'on'));

% Unpack the solution into the exosystem matrices
T = [oG', oq'; qG', qSol'];
[~, o] = sort(T(1, :));
T = T(:, o);
S = [T(2, 1), T(2, 2); T(2, 3), T(2, 4)];
s = [T(2, 5); T(2, 6)];
end

function [S, s] = computeExosystembyFBI2(A, B, x0, q0)

% Enforces S to be full rank
% However, without guessing q0 properly, the eqb point may still be inside
% the simplex

% Solve the system of non-linear equations/constraints
[qSol, ~] = fsolve(@(q) funcFBIConstraints2(q, A, B, x0), ...
    q0, optimset('Display', 'on'));

% Unpack the solution into the exosystem matrices
S = [qSol(1), qSol(2); qSol(3), qSol(4)];
s = [qSol(5); qSol(6)];
end

function bPass = evalIneq(x1, x2, terms, sign, bPrint)

bPass = 0;
dp = dot(x1, x2);
if ~isempty(symvar(dp))
    if bPrint
        fprintf('%s = %s %s 0\n', terms, char(vpa(dp, 4)), sign);
    end
else
    if bPrint
        fprintf('%s = %f %s 0\n', terms, dp, sign);
    end
    
    switch(sign)
        case '<'
            bPass = dp < 0;
        case '<='
            bPass = dp <= 0;
        otherwise
            error('Please use "<" or "<=" only');
    end
end
end

function fNorm = FBISolvedNorm(A, B, C, S, s, Cr, PI, GAMMA, p, gm)

FBI1a = A*PI + B*GAMMA - PI*S;
FBI1b = A*p + B*gm - PI*s;
FBI1c = C*PI - Cr;
FBI1d = C*p;

fNorm = max([norm(FBI1a), norm(FBI1b), norm(FBI1c), norm(FBI1d)]);
end

function [PI, GAMMA, p, gm] = findController(S, s, A, B)

PI = [1 0; S(1,1) S(1,2); 0 1; S(2,1) S(2,2)];
GAMMA = 1/B(2) * [S(1,1)^2 + S(1,2)*S(2,1) - A(2,4)*S(2,1), ...
    S(1,1)*S(1,2) + S(1,2)*S(2,2) - A(2,3) - A(2,4)*S(2,2)];
p = [0; s(1); 0; s(2)];
gm = 1/B(2) * (S(1,1)*s(1) + s(2)*(S(1,2) - A(2,4)));
end

function [xdot] = regulatorODE(t, x, Asys, asys)

xdot = Asys*x + asys;
end