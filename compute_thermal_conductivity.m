% Define a function to compute the thermal conductivity at a given temperature
function conductivity = compute_thermal_conductivity(temperature, coefficients)
    conductivity = coefficients(1) + coefficients(2) * temperature + coefficients(3) * temperature.^2;
end