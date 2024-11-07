% This script utilizes helper functions and resources referenced from Bruinlearn class materials.
% Link: https://bruinlearn.ucla.edu/courses/193842/modules/items/7007232

global Fg M dt
global kappaBar EI GJ voronoiLength
global EA refLen

%% Input parameters

nv = 50; % Number of vertices or nodes in the rod
ne = nv - 1; % Number of edges connecting nodes
ndof = 3 * nv + ne; % Total Degrees of Freedom (DOF)
dt = 0.01; % Step size for time in seconds
totalTime = 5; % Total duration of simulation in seconds
Nsteps = round(totalTime / dt); % Number of iterations in the simulation

RodLength = 0.2; % Length of the rod in meters
natR = 0.02; % Radius of the helix (circle) in meters
r0 = 0.001; % Radius of each node in meters
Y = 10e6; % Young’s modulus 
nu = 0.5; % Poisson’s ratio 
G = Y / (2 * (1 + nu)); % Shear modulus 
rho = 1000; % Density in kg/m^3
g = [0; 0; -9.81]; % Gravitational vector

% Material stiffness values
EI = Y * pi * r0^4 / 4; % Bending stiffness 
GJ = G * pi * r0^4 / 2; % Shearing stiffness
EA = Y * pi * r0^2; % Axial stretching stiffness

% Tolerance value
tol = EI / RodLength^2 * 1e-6;

%% Mass matrix calculation

totalM = pi * r0^2 * RodLength * rho; % Total mass of the rod
dm = totalM / ne; % Uniform mass per edge
massVector = zeros(ndof, 1);

% Assign mass to each node
for c = 1:nv
    ind = [4 * c - 3; 4 * c - 2; 4 * c - 1]; % Node indices
    if c == 1 || c == nv
        massVector(ind) = dm / 2; % Half mass for boundary nodes
    else
        massVector(ind) = dm; % Full mass for interior nodes
    end
end

% Assign rotational inertia to edges
for c = 1:ne
    ind = 4 * c;
    massVector(ind) = 1 / 2 * dm * r0^2;
end

M = diag(massVector); % Create the mass matrix

%% Initialize the DOF vector

nodes = zeros(nv, 3); % Node positions
dTheta = (RodLength / natR) * (1 / ne); % Incremental angle per segment

% Define initial positions based on helical structure
for c = 1:nv
    nodes(c, 1) = natR * cos((c - 1) * dTheta);
    nodes(c, 2) = natR * sin((c - 1) * dTheta);
    nodes(c, 3) = 0;
end

q0 = zeros(ndof, 1); % Initial state vector

% Assign initial positions to the DOF vector
for c = 1:nv
    ind = [4 * c - 3; 4 * c - 2; 4 * c - 1];
    q0(ind) = nodes(c, :);
end

u = zeros(ndof, 1); % Initial velocity vector

%% Compute reference lengths for edges (used for axial stretching)

refLen = zeros(ne, 1);
for c = 1:ne
    dx = nodes(c + 1, :) - nodes(c, :); % Distance between adjacent nodes
    refLen(c) = norm(dx); % Reference length
end

%% Compute Voronoi lengths (associated with each node, for bending and twisting)

voronoiLength = zeros(nv, 1);
for c = 1:nv
    if c == 1
        voronoiLength(c) = 1 / 2 * refLen(c);
    elseif c == nv
        voronoiLength(c) = 1 / 2 * refLen(c - 1);
    else
        voronoiLength(c) = 1 / 2 * refLen(c - 1) + 1 / 2 * refLen(c);
    end
end

%% Reference frame initialization

a1 = zeros(ne, 3); % First reference director
a2 = zeros(ne, 3); % Second reference director

% Compute tangents
tangent = computeTangent(q0);

% Initialize the first frame
t0 = tangent(1, :); % Tangent of the first edge
t1 = [0; 0; -1]; % Arbitrary reference vector
a1Tmp = cross(t0, t1); % Perpendicular vector
if abs(a1Tmp) < 1e-6
    t1 = [0; 1; 0];
    a1Tmp = cross(t0, t1);
end
a1(1, :) = a1Tmp / norm(a1Tmp);
a2(1, :) = cross(tangent(1, :), a1(1, :));

% Propagate reference frame for other edges
for c = 2:ne
    t0 = tangent(c - 1, :);
    t1 = tangent(c, :);
    a1_0 = a1(c - 1, :);
    a1_1 = parallel_transport(a1_0, t0, t1);
    a1(c, :) = a1_1 / norm(a1_1);
    a2(c, :) = cross(t1, a1(c, :));
end

%% Material frame and curvature

theta = q0(4:4:end); % Initial twist angles
[m1, m2] = computeMaterialDirectors(a1, a2, theta);

kappaBar = getkappa(q0, m1, m2); % Initial natural curvature

%% Gravity computation

Fg = zeros(ndof, 1); % Initialize gravity force vector
for c = 1:nv
    ind = [4 * c - 3; 4 * c - 2; 4 * c - 1];
    Fg(ind) = massVector(ind) .* g;
end

%% Time-stepping loop

ctime = 0; % Start time
endZ = zeros(Nsteps, 1); % Record z-coordinate of the last node

for timeStep = 1:Nsteps
    fprintf('Current time = %f\n', ctime);
    ctime = ctime + dt;

    [q, u, a1, a2] = objfun(q0, u, a1, a2, 8:ndof, tol, refTwist);
    q0 = q; % Update DOF vector
    endZ(timeStep) = q(end); % Track end node z-coordinate

    % Plot results at every 5th step
    if mod(timeStep - 1, 5) == 0
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(a1, a2, theta);
        plotrod(q, a1, a2, m1, m2, ctime);
    end
end

% Final displacement plot
figure(2);
timeArray = (1:Nsteps) * dt;
plot(timeArray, endZ, 'b*-'); % Changed to blue with star markers
xlabel('Time (s)');
ylabel('Z-coordinate of last node (m)');
