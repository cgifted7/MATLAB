% Define a function to compute the temperature data based on the current guess for
% the thermal conductivity coefficients
function temperature_data = compute_temperature_data(coefficients)
    data = readmatrix('TwinTech\Day 69\TempDataTest.csv');
    temperature_data = data(1,:); % use initial temperatures as initial condition
    num_times = size(data, 1);
    time_step = 1;
    num_points = size(data, 2);
    density = 1380;
    specific_heat = 983;

    % Compute the temperature at each time step for each point in the material
    for i = 2:num_times
        time = (i - 1) * time_step; % time in seconds
        
        for j = 2:num_points-1
            % Compute the temperature at the current point based on the
            % temperatures at the adjacent points and the thermal conductivity
            % at the current temperature
            temperature = temperature_data(i-1,j);
            conductivity = compute_thermal_conductivity(temperature, coefficients);
            
            left_temp = temperature_data(i-1,j-1);
            right_temp = temperature_data(i-1,j+1);
            heat_flux = conductivity * (right_temp - left_temp) / (2 * time_step);
            heat_transfer = heat_flux / (density * specific_heat);
            
            temperature_data(i,j) = temperature + heat_transfer;
        end
    end
end