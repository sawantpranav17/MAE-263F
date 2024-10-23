% Euler Bernoulli Result: 0.1122937432 m
% Portion of this code and the helper functions/files have been referenced
% From the modules posted on Bruinlearn
% Link: https://colab.research.google.com/drive/1mv5eBdCb42CvJXPfXmdaKfEIJdb5fgNw?usp=sharing
N_nodes = 50; % Number of Nodes

num_dof = N_nodes * 2; % Degrees of Freedom
time_step = 0.001; % Time step

rod_length = 1; % Total length of the rod
delta_length = rod_length / (N_nodes - 1); % Distance between each node

load_force = 2000; % Load in Newtons
distance_load = 0.75; % Distance at which the load is acting
load_node = round(distance_load / delta_length);

outer_radius = 0.013; % Outer Radius
inner_radius = 0.011; % Inner Radius

% Material Properties
density_aluminum = 2700; % Density of Spheres
gravity = 9.8; % Acceleration due to gravity

total_simulation_time = 1; % Total time of simulation

% Utility Parameters
num_edges = N_nodes - 1; % Number of Edges

youngs_modulus = 70 * 1e9;
moment_of_inertia = pi * (outer_radius^4 - inner_radius^4) / 4;
cross_section_area = pi * (outer_radius^2 - inner_radius^2);
bending_stiffness = youngs_modulus * moment_of_inertia;
tensile_stiffness = youngs_modulus * cross_section_area;

% Initial Configuration
node_positions = zeros(N_nodes, 2);
for i = 1:N_nodes
    node_positions(i, 1) = (i - 1) * delta_length; % x-coordinates
    node_positions(i, 2) = 0; % Can be removed
end

% Mass Matrix
mass_matrix = zeros(num_dof, num_dof);
mass_i = pi * (outer_radius^2 - inner_radius^2) * rod_length * density_aluminum / (N_nodes - 1);
for i = 1:N_nodes
    mass_matrix(2 * i - 1, 2 * i - 1) = mass_i;
    mass_matrix(2 * i, 2 * i) = mass_i;
end

% Load Matrix
load_vector = zeros(num_dof, 1);
load_vector(2 * load_node) = load_vector(2 * load_node) - load_force;

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

max_deflection = zeros(num_steps, 1); % Maximum deflection of the beam

for step = 2:num_steps
    fprintf('Time = %f\n', (step - 1) * time_step);

    fixed_dof = [1; 2; num_dof];
    free_dof = 3:num_dof - 1;
    boundary_condition_vector = [0; 0; 0];

    % Guess
    current_positions = initial_positions; % Setting the new positions to the old ones
    current_positions(fixed_dof) = boundary_condition_vector;

    % Newton-Raphson Method
    error = 10 * tolerance;
    while error > tolerance
        free_positions = current_positions(free_dof);

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

        % Load
        force_vector = force_vector - load_vector;

        free_force_vector = force_vector(free_dof);
        free_jacobian_matrix = jacobian_matrix(free_dof, free_dof);

        dq_free = free_jacobian_matrix \ free_force_vector;

        % Update
        free_positions = free_positions - dq_free;
        error = sum(abs(free_force_vector));

        current_positions(free_dof) = free_positions;
    end

    % New Velocity
    velocity_vector = (current_positions - initial_positions) / time_step;
    initial_positions = current_positions;

    % Calculating Max Deflection
    max_deflection_value = 0;
    for i = 1:N_nodes
        if max_deflection_value > current_positions(2 * i)
            max_deflection_value = current_positions(2 * i);
        end
    end

    % Storing Max Deflection
    max_deflection(step) = max_deflection_value;

    % Plot
    figure(1);
    plot(current_positions(1:2:end), current_positions(2:2:end), 'ro-');
    axis equal
    xlabel('x (meter)')
    ylabel('y (meter)')
    drawnow
end

figure(2);
time_array = (1:num_steps) * time_step;
plot(time_array, max_deflection, 'k-');
xlabel('Time, t [sec]');
ylabel('Maximum Deflection, y [meter]');


% Code for Question 02

% %Portion of this code and the helper functions/files have been referenced
% %from the class resources on Bruinlearn
% %Link: https://bruinlearn.ucla.edu/courses/171817/modules
% 
% load_array = zeros(25);
% max_deflection_array = zeros(25);
% for loop = 1:25
%     % Euler Bernoulli Result: 0.1122937432 m
% 
%     N_nodes = 50; % Number of Nodes
%     num_dof = N_nodes * 2; % Degrees of Freedom
%     time_step = 0.01; % Time step
% 
%     rod_length = 1; % Total length of the rod
% delta_length = rod_length / (N_nodes - 1); % Distance between each node
% 
%     load_force = 1000 + (loop - 1) * 1000; % Load in Newtons
%     fprintf('Load = %f\n', load_force);
% 
%     load_array(loop) = load_force;
%     distance_load = 0.75; % Distance at which the load is acting
%     load_node = round(distance_load / delta_length);
% 
%     outer_radius = 0.013; % Outer Radius
%     inner_radius = 0.011; % Inner Radius
% 
%     % Material Properties
%     density_aluminum = 2700; % Density of Spheres
%     gravity = 9.8; % Acceleration due to gravity
% 
%     total_simulation_time = 1; % Total time of simulation
% 
%     % Utility Parameters
%     num_edges = N_nodes - 1; % Number of Edges
% 
% youngs_modulus = 70 * 1e9;
% moment_of_inertia = pi * (outer_radius^4 - inner_radius^4) / 4;
% cross_section_area = pi * (outer_radius^2 - inner_radius^2);
% bending_stiffness = youngs_modulus * moment_of_inertia;
% tensile_stiffness = youngs_modulus * cross_section_area;
% 
%     % Initial Configuration
%     node_positions = zeros(N_nodes, 2);
%     for i = 1:N_nodes
%         node_positions(i, 1) = (i - 1) * delta_length; % x-coordinates
%         node_positions(i, 2) = 0; % Can be removed
%     end
% 
%     % Mass Matrix
%     mass_matrix = zeros(num_dof, num_dof);
%     mass_i = pi * (outer_radius^2 - inner_radius^2) * rod_length * density_aluminum / (N_nodes - 1);
%     for i = 1:N_nodes
%         mass_matrix(2 * i - 1, 2 * i - 1) = mass_i;
%         mass_matrix(2 * i, 2 * i) = mass_i;
%     end
% 
%     % Load Matrix
%     load_vector = zeros(num_dof, 1);
%     load_vector(2 * load_node) = load_vector(2 * load_node) - load_force;
% 
%     % Initial Positions
%     initial_positions = zeros(num_dof, 1);
%     for i = 1:N_nodes
%         initial_positions(2 * i - 1) = node_positions(i, 1);
%         initial_positions(2 * i) = node_positions(i, 2);
%     end
% 
%     % Initial Velocities
%     current_positions = initial_positions; % Degrees of Freedom vector
%     velocity_vector = (current_positions - initial_positions) / time_step; % Velocity vector
% 
%     % Tolerance
%     tolerance = bending_stiffness / rod_length^2 * 1e-3;
% 
%     num_steps = round(total_simulation_time / time_step);
% 
%     max_deflection = zeros(num_steps, 1); % Maximum deflection of the beam
% 
%     for step = 2:num_steps
%         %fprintf('Time = %f\n', (step - 1) * time_step);
% 
%         fixed_dof = [1; 2; num_dof];
%         free_dof = 3:num_dof - 1;
%         boundary_condition_vector = [0; 0; 0];
% 
%         % Guess
%         current_positions = initial_positions; % Setting the new positions to the old ones
%         current_positions(fixed_dof) = boundary_condition_vector;
% 
%         % Newton-Raphson Method
%         error = 10 * tolerance;
%         while error > tolerance
%             free_positions = current_positions(free_dof);
% 
%             force_vector = mass_matrix / time_step * ((current_positions - initial_positions) / time_step - velocity_vector);
%             jacobian_matrix = mass_matrix / time_step^2;
% 
%             % Effects of forces
%             for i = 1:N_nodes - 1
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
%             for i = 2:N_nodes - 1
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
%             % Load
%             force_vector = force_vector - load_vector;
% 
%             free_force_vector = force_vector(free_dof);
%             free_jacobian_matrix = jacobian_matrix(free_dof, free_dof);
% 
%             dq_free = free_jacobian_matrix \ free_force_vector;
% 
%             % Update
%             free_positions = free_positions - dq_free;
%             error = sum(abs(free_force_vector));
% 
%             current_positions(free_dof) = free_positions;
%         end
% 
%         % New Velocity
%         velocity_vector = (current_positions - initial_positions) / time_step;
%         initial_positions = current_positions;
% 
%         % Calculating Max Deflection
%         max_deflection_value = 0;
%         for i = 1:N_nodes
%             if max_deflection_value > current_positions(2 * i)
%                 max_deflection_value = current_positions(2 * i);
%             end
%         end
% 
%         % Storing Max Deflection
%         max_deflection(step) = max_deflection_value;
%     end
%     max_deflection_array(loop) = max_deflection(num_steps - 1);
% end
% 
% % Euler-Bernoulli Deflection Calculation
% euler_bernoulli_deflection = zeros(25, 1);
% for loop = 1:25
%     euler_bernoulli_deflection(loop) = -load_array(loop) * 0.25 * (1 - 0.25 * 0.25)^(1.5) / (9 * sqrt(3) * youngs_modulus * moment_of_inertia);
% end
% 
% % Plotting Results
% figure(1);
% plot(load_array, max_deflection_array);
% hold on;
% plot(load_array, euler_bernoulli_deflection, '-r');
% xlabel('Load Applied, Newtons');
% ylabel('Maximum Deflection, meters');
% legend('Simulation', 'Euler–Bernoulli beam theory');
% hold off;
% 
