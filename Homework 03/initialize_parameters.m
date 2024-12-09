%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252
%%File Name: initialize_parameters.m

function parameters = initialize_parameters(layer_dims)
% Initializes the weights and biases of the feedforward neural network
% Inputs:
%   layer_dims: array of layer dimensions, including input and output layers
% Output: 
%   parameters: a struct containing W1, b1, W2, b2, etc.
%
% Hint: parameters{1}.W should be a ( layer_dims(2) x 784 ) matrix, which 
% is the initialized weights connecting the input layer and the first 
% hidden layer. parameters{1}.b should be a ( layer_dims(2) x 1 ) array, 
% which is the initialized biases connecting the input layer and the first 
% hidden layer.

    L = length(layer_dims);
    parameters = cell(1, L-1);
    for l = 1:(L-1)
        parameters{l}.W  = randn(layer_dims(l+1), layer_dims(l));
        parameters{l}.b = zeros(layer_dims(l+1), 1);
    end
end