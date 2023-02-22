% Load the temperature data from the CSV file
data = readmatrix('TempDataTest');

% Set up the problem parameters
density = 1380; % kg/m^3
specific_heat = 983; % J/(kg*K)
time_step = 1; % seconds
num_points = size(data, 2); % number of temperature data points
num_times = size(data, 1); % number of time steps
initial_guess = [1.187*10^-6, -.0012649, 0.87]; % initial guess for thermal conductivity coefficients
tolerance = 1e-5; % tolerance for convergence
max_iterations = 1000; % maximum number of iterations

% Implement the Levenberg-Marquardt algorithm to solve for the thermal
% conductivity coefficients
coefficients = initial_guess;
lambda = 0.01;

for i = 1:max_iterations
    % Compute the Jacobian matrix
    J = zeros(num_times * num_points, 3);
    for j = 1:num_points-1
        temperature = data(:,j);
        conductivity = compute_thermal_conductivity(temperature, coefficients);
        
        for k = 2:num_times-1
            left_temp = temperature(k-1);
            
            right_temp = temperature(k+1);
          
            heat_flux = conductivity(k) * (right_temp - left_temp) / (2 * time_step);
            J((k-1)*num_points+j,:) = [-heat_flux/temperature(k), -heat_flux*temperature(k)/temperature(k+1), -heat_flux*temperature(k)*temperature(k)/temperature(k+1)^2];
        end
    end
    
    % Compute the current error and gradient
    error = objective_function(coefficients);
    gradient = J' * error;
    
    % Compute the new guess for the coefficients
    new_coefficients = coefficients - (J' * J + lambda * eye(3)) \ (J' * error);
    disp(new_coefficients);

    % Compute the new error and check for convergence
    new_error = objective_function(new_coefficients);
    if norm(new_error) < tolerance && norm(new_coefficients - coefficients) < tolerance
        coefficients = new_coefficients;
        break;
    end

    % Update the guess and the damping parameter
    if norm(new_error) < norm(error)
        coefficients = new_coefficients;
        lambda = lambda / 10;
    else
        lambda = lambda * 10;
    end
end

% Print the final guess for the thermal conductivity coefficients
fprintf('Final guess for the thermal conductivity coefficients: %f, %f, %f\n', coefficients(1), coefficients(2), coefficients(3));

% Compute the predicted temperature data using the final guess for the thermal
% conductivity coefficients
predicted_data = compute_temperature_data(coefficients);

% Plot the observed and predicted temperature data
figure;
plot(0:time_step:(num_times-1)*time_step, data, '-');
%hold on;
figure;
plot(0:time_step:(num_times-1)*time_step, predicted_data, '--');
xlabel('Time (s)');
ylabel('Temperature (C)');
legend('Observed', 'Predicted');


