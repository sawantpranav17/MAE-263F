% Portion of this code and the helper functions/files have been referenced
% From the modules posted on Bruinlearn
% Link: https://bruinlearn.ucla.edu/courses/193842/files/18268548?module_item_id=6975393

N_nodes = 21; % Number of Nodes

num_dof = N_nodes * 2; % Degrees of Freedom
time_step = 0.01; % Time step

rod_length = 0.1; % Total length of the rod
delta_length = rod_length / (N_nodes - 1); % Distance between each node

% Radii of Spheres
radius = zeros(N_nodes, 1);
mid_node = (N_nodes + 1) / 2;
radius(:) = delta_length / 10;
radius(mid_node) = 0.025;

% Material Properties
density_spheres = 7000; % Density of Spheres
density_fluid = 1000; % Density of Fluid
delta_density = density_spheres - density_fluid; % Difference between densities
viscosity = 1000; % Viscosity

rod_radius = 0.001; % Radius of the rod
youngs_modulus = 1e9; % Young's Modulus
gravity = 9.8; % Acceleration due to gravity

total_simulation_time = 50; % Total time of simulation

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
for i = 1:N_nodes
    mass_i = 4 / 3 * pi * radius(i)^3 * density_spheres;
    mass_matrix(2 * i - 1, 2 * i - 1) = mass_i;
    mass_matrix(2 * i, 2 * i) = mass_i;
end

% Damping Matrix
damping_matrix = zeros(num_dof, num_dof);
for i = 1:N_nodes
    damping_i = 6 * pi * viscosity * radius(i);
    damping_matrix(2 * i - 1, 2 * i - 1) = damping_i;
    damping_matrix(2 * i, 2 * i) = damping_i;
end

% Weight Matrix
weight_vector = zeros(num_dof, 1);
for i = 1:N_nodes
    weight_vector(2 * i) = -4 / 3 * pi * radius(i)^3 * delta_density * gravity;
end

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
mid_node_y(1) = current_positions(2 * mid_node);
mid_node_v(1) = velocity_vector(2 * mid_node);

for step = 2:num_steps
    fprintf('Time = %f\n', (step - 1) * time_step);

    % Guess
    current_positions = initial_positions; % Setting the new positions to the old ones
    % Newton-Raphson Method
    error = 10 * tolerance;
    while error > tolerance
        force_vector = mass_matrix / time_step * ((current_positions - initial_positions) / time_step - velocity_vector);
        jacobian_matrix = mass_matrix / time_step^2;

        % Effects of forces
        for i = 1:N_nodes - 1
            x_i = current_positions(2 * i - 1);
            y_i = current_positions(2 * i);
            x_ip1 = current_positions(2 * i + 1);
            y_ip1 = current_positions(2 * i + 2);
            segment_length = delta_length;
            force_increment = gradEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
            jacobian_increment = hessEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
            index = [2 * i - 1, 2 * i, 2 * i + 1, 2 * i + 2];
            force_vector(index) = force_vector(index) + force_increment;
            jacobian_matrix(index, index) = jacobian_matrix(index, index) + jacobian_increment;
        end

        for i = 2:N_nodes - 1
            x_im1 = current_positions(2 * i - 3);
            y_im1 = current_positions(2 * i - 2);
            x_i = current_positions(2 * i - 1);
            y_i = current_positions(2 * i);
            x_ip1 = current_positions(2 * i + 1);
            y_ip1 = current_positions(2 * i + 2);
            curvature = 0;
            segment_length = delta_length;
            force_increment = gradEb(x_im1, y_im1, x_i, y_i, x_ip1, y_ip1, curvature, segment_length, bending_stiffness);
            jacobian_increment = hessEb(x_im1, y_im1, x_i, y_i, x_ip1, y_ip1, curvature, segment_length, bending_stiffness);
            index = [2 * i - 3, 2 * i - 2, 2 * i - 1, 2 * i, 2 * i + 1, 2 * i + 2];
            force_vector(index) = force_vector(index) + force_increment;
            jacobian_matrix(index, index) = jacobian_matrix(index, index) + jacobian_increment;
        end


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

    mid_node_y(step) = current_positions(2 * mid_node);
    mid_node_v(step) = velocity_vector(2 * mid_node);
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

mid_node_v(num_steps - 1)




% code for question 3

% node_count_array = zeros(19);
% terminal_velocity_array = zeros(19);
% 
% time_step = 0.01; % Time step
% 
% rod_length = 0.1; % Total length of the rod
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
% total_simulation_time = 50; % Total time of simulation
% 
% bending_stiffness = youngs_modulus * pi * rod_radius^4 / 4; % Bending Stiffness
% tensile_stiffness = youngs_modulus * pi * rod_radius^2;
% 
% % Tolerance
% tolerance = bending_stiffness / rod_length^2 * 1e-3;
% 
% num_steps = round(total_simulation_time / time_step);
% 
% all_mid_y = zeros(num_steps, 1); % y-position of middle node
% all_mid_v = zeros(num_steps, 1); % y-velocity of middle node
% 
% for loop = 1:19
%     node_count = 3 + (loop - 1) * 2;
%     fprintf('Number of Nodes = %f\n', node_count);
%     node_count_array(loop) = node_count;
%     num_dof = node_count * 2; % Degrees of Freedom
% 
%     delta_length = rod_length / (node_count - 1); % Distance between each node
% 
%     % Radii of Spheres
%     radius = zeros(node_count, 1);
%     mid_node = (node_count + 1) / 2;
%     radius(:) = delta_length / 10;
%     radius(mid_node) = 0.025;
% 
%     % Utility Parameters
%     num_edges = node_count - 1; % Number of Edges
% 
%     % Initial Configuration
%     node_positions = zeros(node_count, 2);
%     for i = 1:node_count
%         node_positions(i, 1) = (i - 1) * delta_length; % x-coordinates
%         node_positions(i, 2) = 0; % Can be removed
%     end
% 
%     % Mass Matrix
%     mass_matrix = zeros(num_dof, num_dof);
%     for i = 1:node_count
%         mass_i = 4 / 3 * pi * radius(i)^3 * density_spheres;
%         mass_matrix(2 * i - 1, 2 * i - 1) = mass_i;
%         mass_matrix(2 * i, 2 * i) = mass_i;
%     end
% 
%     % Damping Matrix
%     damping_matrix = zeros(num_dof, num_dof);
%     for i = 1:node_count
%         damping_i = 6 * pi * viscosity * radius(i);
%         damping_matrix(2 * i - 1, 2 * i - 1) = damping_i;
%         damping_matrix(2 * i, 2 * i) = damping_i;
%     end
% 
%     % Weight Matrix
%     weight_vector = zeros(num_dof, 1);
%     for i = 1:node_count
%         weight_vector(2 * i) = -4 / 3 * pi * radius(i)^3 * delta_density * gravity;
%     end
% 
%     % Initial Positions
%     initial_positions = zeros(num_dof, 1);
%     for i = 1:node_count
%         initial_positions(2 * i - 1) = node_positions(i, 1);
%         initial_positions(2 * i) = node_positions(i, 2);
%     end
% 
%     % Initial Velocities
%     current_positions = initial_positions; % Degrees of Freedom vector
%     velocity_vector = (current_positions - initial_positions) / time_step; % Velocity vector
% 
%     for step = 2:num_steps
%         %fprintf('Time = %f\n', (step - 1) * time_step);
% 
%         % Guess
%         current_positions = initial_positions; % Setting the new positions to the old ones
%         % Newton-Raphson Method
%         error = 10 * tolerance;
%         while error > tolerance
%             force_vector = mass_matrix / time_step * ((current_positions - initial_positions) / time_step - velocity_vector);
%             jacobian_matrix = mass_matrix / time_step^2;
% 
%             % Effects of forces
%             for i = 1:node_count - 1
%                 x_i = current_positions(2 * i - 1);
%                 y_i = current_positions(2 * i);
%                 x_ip1 = current_positions(2 * i + 1);
%                 y_ip1 = current_positions(2 * i + 2);
%                 segment_length = delta_length;
%                 force_increment = gradEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
%                 jacobian_increment = hessEs(x_i, y_i, x_ip1, y_ip1, segment_length, tensile_stiffness);
%                 index = [2 * i - 1, 2 * i, 2 * i + 1, 2 * i + 2];
%                 force_vector(index) = force_vector(index) + force_increment;
%                 jacobian_matrix(index, index) = jacobian_matrix(index, index) + jacobian_increment;
%             end
% 
%             for i = 2:node_count - 1
%                 x_im1 = current_positions(2 * i - 3);
%                 y_im1 = current_positions(2 * i - 2);
%                 x_i = current_positions(2 * i - 1);
%                 y_i = current_positions(2 * i);
%                 x_ip1 = current_positions(2 * i + 1);
%                 y_ip1 = current_positions(2 * i + 2);
%                 curvature = 0;
%                 segment_length = delta_length;
%                 force_increment = gradEb(x_im1, y_im1, x_i, y_i, x_ip1, y_ip1, curvature, segment_length, bending_stiffness);
%                 jacobian_increment = hessEb(x_im1, y_im1, x_i, y_i, x_ip1, y_ip1, curvature, segment_length, bending_stiffness);
%                 index = [2 * i - 3, 2 * i - 2, 2 * i - 1, 2 * i, 2 * i + 1, 2 * i + 2];
%                 force_vector(index) = force_vector(index) + force_increment;
%                 jacobian_matrix(index, index) = jacobian_matrix(index, index) + jacobian_increment;
%             end
% 
%             % Viscous Force
%             force_vector = force_vector + damping_matrix * (current_positions - initial_positions) / time_step;
%             jacobian_matrix = jacobian_matrix + damping_matrix / time_step;
% 
%             % Weight
%             force_vector = force_vector - weight_vector;
% 
%             % Update
%             current_positions = current_positions - jacobian_matrix \ force_vector;
%             error = sum(abs(force_vector));
%         end
% 
%         % New Velocity
%         velocity_vector = (current_positions - initial_positions) / time_step;
%         initial_positions = current_positions;
% 
%         all_mid_y(step) = current_positions(2 * mid_node);
%         all_mid_v(step) = velocity_vector(2 * mid_node);
%     end
% 
%     terminal_velocity_array(loop) = all_mid_v(num_steps - 1);
% end
% 
% figure(1);
% plot(node_count_array, terminal_velocity_array, 'k-');
% axis([0 45 -0.006 -0.0058]);
% xlabel('Number of Spheres, N');
% ylabel('Terminal Velocity, v [meter/sec]');
% 
