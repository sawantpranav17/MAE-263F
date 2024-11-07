% Portion of this code and the helper functions/files have been referenced
% From the modules posted on Bruinlearn
% https://bruinlearn.ucla.edu/courses/193842/files/18268541?module_item_id=6975388


N_nodes = 3; % Number of Nodes

num_dof = N_nodes * 2; % Degrees of Freedom
time_step = 0.01; % Time step

rod_length = 0.1; % Total length of the rod
delta_length = rod_length / (N_nodes - 1); % Distance between each node

% Radii of Spheres
radius1 = 0.005;
radius2 = 0.025;
radius3 = 0.005;

% Material Properties
density_spheres = 7000; % Density of Spheres
density_fluid = 1000; % Density of Fluid
delta_density = density_spheres - density_fluid; % Difference between densities
viscosity = 1000; % Viscosity

rod_radius = 0.001; % Radius of the rod
youngs_modulus = 1e9; % Young's Modulus
gravity = 9.8; % Acceleration due to gravity

total_simulation_time = 10; % Total time of simulation

% Utility Parameters
num_edges = N_nodes - 1; % Number of Edges
bending_stiffness = youngs_modulus * pi * rod_radius^4 / 4; % Bending Stiffness
tensile_stiffness = youngs_modulus * pi * rod_radius^2;

% Initial Configuration
node_positions = zeros(N_nodes, 2);
for i = 1:N_nodes
    node_positions(i, 1) = (i - 1) * delta_length; % x-coordinates
    node_positions(i, 2) = 0; % Can be removed
end

% Mass Matrix
mass_matrix = zeros(num_dof, num_dof);
mass_matrix(1, 1) = 4 / 3 * pi * radius1^3 * density_spheres;
mass_matrix(2, 2) = 4 / 3 * pi * radius1^3 * density_spheres;
mass_matrix(3, 3) = 4 / 3 * pi * radius2^3 * density_spheres;
mass_matrix(4, 4) = 4 / 3 * pi * radius2^3 * density_spheres;
mass_matrix(5, 5) = 4 / 3 * pi * radius3^3 * density_spheres;
mass_matrix(6, 6) = 4 / 3 * pi * radius3^3 * density_spheres;

% Damping Matrix
damping_matrix = zeros(num_dof, num_dof);
damping_matrix(1, 1) = 6 * pi * viscosity * radius1;
damping_matrix(2, 2) = 6 * pi * viscosity * radius1;
damping_matrix(3, 3) = 6 * pi * viscosity * radius2;
damping_matrix(4, 4) = 6 * pi * viscosity * radius2;
damping_matrix(5, 5) = 6 * pi * viscosity * radius3;
damping_matrix(6, 6) = 6 * pi * viscosity * radius3;

% Weight Matrix
weight_vector = zeros(num_dof, 1);
weight_vector(2) = -4 / 3 * pi * radius1^3 * delta_density * gravity;
weight_vector(4) = -4 / 3 * pi * radius2^3 * delta_density * gravity;
weight_vector(6) = -4 / 3 * pi * radius3^3 * delta_density * gravity;

% Initial Positions
initial_positions = zeros(num_dof, 1);
for i = 1:N_nodes
    initial_positions(2 * i - 1) = node_positions(i, 1);
    initial_positions(2 * i) = node_positions(i, 2);
end

% Initial Velocities
current_positions = initial_positions; % Degrees of Freedom vector
velocity_vector = (current_positions - initial_positions) / time_step; % Velocity vector

% Tolerance
tolerance = bending_stiffness / rod_length^2 * 1e-3;

num_steps = round(total_simulation_time / time_step);

mid_node_y = zeros(num_steps, 1); % y-position of middle node
mid_node_v = zeros(num_steps, 1); % y-velocity of middle node
mid_node_y(1) = current_positions(4);
mid_node_v(1) = velocity_vector(4);

% Times of interest
times_of_interest = [0, 0.01, 0.05, 0.1, 1.0, 10.0];
positions_of_interest = [];

for step = 2:num_steps
    current_time = (step - 1) * time_step;
    fprintf('Time = %f\n', current_time);

    % Guess
    current_positions = initial_positions; % Setting the new positions to the old ones
    % Newton-Raphson Method
    error = 10 * tolerance;
    while error > tolerance
        force_vector = mass_matrix / time_step * ((current_positions - initial_positions) / time_step - velocity_vector);
        jacobian_matrix = mass_matrix / time_step^2;

        % Effects of forces

        % Between Nodes 1 and 2
        x_i = current_positions(1);
        y_i = current_positions(2);
        x_ip1 = current_positions(3);
        y_ip1 = current_positions(4);
        segment_length = delta_length;
        force_increment = gradEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
        jacobian_increment = hessEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
        force_vector(1:4) = force_vector(1:4) + force_increment;
        jacobian_matrix(1:4, 1:4) = jacobian_matrix(1:4, 1:4) + jacobian_increment;

        % Between Nodes 2 and 3
        x_i = current_positions(3);
        y_i = current_positions(4);
        x_ip1 = current_positions(5);
        y_ip1 = current_positions(6);
        segment_length = delta_length;
        force_increment = gradEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
        jacobian_increment = hessEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
        force_vector(3:6) = force_vector(3:6) + force_increment;
        jacobian_matrix(3:6, 3:6) = jacobian_matrix(3:6, 3:6) + jacobian_increment;

        % Bending at Node 2
        x_im1 = current_positions(1);
        y_im1 = current_positions(2);
        x_i = current_positions(3);
        y_i = current_positions(4);
        x_ip1 = current_positions(5);
        y_ip1 = current_positions(6);
        curvature = 0;
        segment_length = delta_length;
        force_increment = gradEb(x_im1, y_im1, x_i, y_i, x_ip1, y_ip1, curvature, segment_length, bending_stiffness);
        jacobian_increment = hessEb(x_im1, y_im1, x_i, y_i, x_ip1, y_ip1, curvature, segment_length, bending_stiffness);
        force_vector(1:6) = force_vector(1:6) + force_increment;
        jacobian_matrix(1:6, 1:6) = jacobian_matrix(1:6, 1:6) + jacobian_increment;

        % Viscous Force
        force_vector = force_vector + damping_matrix * (current_positions - initial_positions) / time_step;
        jacobian_matrix = jacobian_matrix + damping_matrix / time_step;

        % Weight
        force_vector = force_vector - weight_vector;

        % Update
        current_positions = current_positions - jacobian_matrix \ force_vector;
        error = sum(abs(force_vector));
    end

    % New Velocity
    velocity_vector = (current_positions - initial_positions) / time_step;
    initial_positions = current_positions;

    % Plot
    figure(1);
    plot(current_positions(1:2:end), current_positions(2:2:end), 'ko-');
    axis equal
    xlabel('x (meter)')
    ylabel('y (meter)')
    drawnow

    mid_node_y(step) = current_positions(4);
    mid_node_v(step) = velocity_vector(4);

    % Store positions at times of interest
    if any(abs(current_time - times_of_interest) < 1e-5)
        positions_of_interest = [positions_of_interest; current_positions'];
    end
end

% Plotting the position
figure(2);
time_array = (1:num_steps) * time_step;
plot(time_array, mid_node_y, 'k-');
xlabel('Time, t [sec]');
ylabel('Position of middle node, y [meter]');

% Plotting the velocity
figure(3);
plot(time_array, mid_node_v, 'k-');
xlabel('Time, t [sec]');
ylabel('Velocity of middle node, v [meter/sec]');

% Plotting positions at specific times
figure(4);
hold on;
colors = ['r', 'b', 'g', 'c', 'm', 'y'];
legends = {};
for i = 1:min(size(positions_of_interest, 1), length(colors))
    plot(positions_of_interest(i, 1:2:end), positions_of_interest(i, 2:2:end), 'o-', 'Color', colors(i));
    legends{i} = ['t = ', num2str(times_of_interest(i))];
end
hold off;
xlabel('x (meter)');
ylabel('y (meter)');
legend(legends);



% Below code is for explicit method

% Portion of this code have been referenced
% From the modules posted on Bruinlearn
% https://bruinlearn.ucla.edu/courses/193842/files/18268541?module_item_id=6975388


% N_nodes = 3; % Number of Nodes
% 
% num_dof = N_nodes * 2; % Degrees of Freedom
% time_step = 0.00001; % Time step
% 
% rod_length = 1; % Total length of the rod
% delta_length = rod_length / (N_nodes - 1); % Distance between each node
% 
% % Radii of Spheres
% radius1 = 0.005;
% radius2 = 0.025;
% radius3 = 0.005;
% 
% % Material Properties
% density_spheres = 7000; % Density of Spheres
% density_fluid = 1000; % Density of Fluid
% delta_density = density_spheres - density_fluid; % Difference between densities
% viscosity = 1000; % Viscosity
% 
% rod_radius = 0.001; % Radius of the rod
% youngs_modulus = 1e9; % Young's Modulus
% gravity = 9.8; % Acceleration due to gravity
% 
% total_simulation_time = 10; % Total time of simulation
% 
% % Utility Parameters
% num_edges = N_nodes - 1; % Number of Edges
% bending_stiffness = youngs_modulus * pi * rod_radius^4 / 4; % Bending Stiffness
% tensile_stiffness = youngs_modulus * pi * rod_radius^2;
% 
% % Initial Configuration
% node_positions = zeros(N_nodes, 2);
% for i = 1:N_nodes
%     node_positions(i, 1) = (i - 1) * delta_length; % x-coordinates
%     node_positions(i, 2) = 0; % Can be removed
% end
% 
% % Mass Matrix
% mass_matrix = zeros(num_dof, num_dof);
% mass_matrix(1, 1) = 4 / 3 * pi * radius1^3 * density_spheres;
% mass_matrix(2, 2) = 4 / 3 * pi * radius1^3 * density_spheres;
% mass_matrix(3, 3) = 4 / 3 * pi * radius2^3 * density_spheres;
% mass_matrix(4, 4) = 4 / 3 * pi * radius2^3 * density_spheres;
% mass_matrix(5, 5) = 4 / 3 * pi * radius3^3 * density_spheres;
% mass_matrix(6, 6) = 4 / 3 * pi * radius3^3 * density_spheres;
% 
% % Damping Matrix
% damping_matrix = zeros(num_dof, num_dof);
% damping_matrix(1, 1) = 6 * pi * viscosity * radius1;
% damping_matrix(2, 2) = 6 * pi * viscosity * radius1;
% damping_matrix(3, 3) = 6 * pi * viscosity * radius2;
% damping_matrix(4, 4) = 6 * pi * viscosity * radius2;
% damping_matrix(5, 5) = 6 * pi * viscosity * radius3;
% damping_matrix(6, 6) = 6 * pi * viscosity * radius3;
% 
% % Weight Matrix
% weight_vector = zeros(num_dof, 1);
% weight_vector(2) = -4 / 3 * pi * radius1^3 * delta_density * gravity;
% weight_vector(4) = -4 / 3 * pi * radius2^3 * delta_density * gravity;
% weight_vector(6) = -4 / 3 * pi * radius3^3 * delta_density * gravity;
% 
% % Initial Positions
% initial_positions = zeros(num_dof, 1);
% for i = 1:N_nodes
%     initial_positions(2 * i - 1) = node_positions(i, 1);
%     initial_positions(2 * i) = node_positions(i, 2);
% end
% 
% % Initial Velocities
% current_positions = initial_positions; % Degrees of Freedom vector
% velocity_vector = (current_positions - initial_positions) / time_step; % Velocity vector
% 
% % Tolerance
% tolerance = bending_stiffness / rod_length^2 * 1e-3;
% 
% num_steps = round(total_simulation_time / time_step);
% 
% mid_node_y = zeros(num_steps, 1); % y-position of middle node
% mid_node_v = zeros(num_steps, 1); % y-velocity of middle node
% mid_node_y(1) = current_positions(4);
% mid_node_v(1) = velocity_vector(4);
% 
% for step = 2:num_steps
%     fprintf('Time = %f\n %f\n', (step - 1) * time_step, current_positions(4));
%     force_vector = -damping_matrix * velocity_vector + weight_vector;
%     force_vector(1:4) = force_vector(1:4) - gradEs(current_positions(1), current_positions(2), current_positions(3), current_positions(4), l_k, tensile_stiffness);
%     force_vector(3:6) = force_vector(3:6) - gradEs(current_positions(3), current_positions(4), current_positions(5), current_positions(6), l_k, tensile_stiffness);
%     force_vector(1:6) = force_vector(1:6) - gradEb(current_positions(1), current_positions(2), current_positions(3), current_positions(4), current_positions(5), current_positions(6), curv, l_k, bending_stiffness);
%     % New velocity
%     for i = 1:N_nodes
%         current_positions(2 * i - 1) = initial_positions(2 * i - 1) + time_step * (velocity_vector(2 * i - 1) + force_vector(2 * i - 1) * time_step / mass_matrix(2 * i - 1, 2 * i - 1));
%         current_positions(2 * i) = initial_positions(2 * i) + time_step * (velocity_vector(2 * i) + force_vector(2 * i) * time_step / mass_matrix(2 * i, 2 * i));
%     end
% 
%     % Plot
%     figure(1);
%     plot(current_positions(1:2:end), current_positions(2:2:end), 'ko-');
%     axis equal
%     xlabel('x (meter)')
%     ylabel('y (meter)')
%     drawnow
% 
%     velocity_vector = (current_positions - initial_positions) / time_step;
%     initial_positions = current_positions;
% end
% figure(1);
% plot(current_positions(1:2:end), current_positions(2:2:end), 'ko-');
% axis equal
% xlabel('x (meter)')
% ylabel('y (meter)')
% drawnow
