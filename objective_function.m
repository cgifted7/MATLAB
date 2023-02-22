% Define the objective function
function error = objective_function(coefficients)
    % Compute the predicted temperature data based on the current guess for the
    % thermal conductivity coefficients
    data = readmatrix('TwinTech\Day 69\TempDataTest.csv');
    predicted_data = compute_temperature_data(coefficients);
    
    % Compute the difference between the observed temperature data and the
    % predicted temperature data
    error = data - predicted_data;
    
    % Flatten the error matrix into a vector
    error = error(:);
end