%This code and the helper functions/files have been referenced
%from the class resources on Bruinlearn

%Name: Pranav Sawant
%%UID: 906146252
%%File Name: backward_propagation.m

function gradients = backward_propagation(X, Y, parameters, activations)
% Computes the gradients of a feedforward neural network given input data, true labels, and learned parameters
% X: input data, shape (input size, number of examples)
% Y: true labels, shape (output size, number of examples)
% parameters: learned parameters, a struct containing W1, b1, W2, b2, etc.
% activations: array of activations at each layer, including input and output layers, computed using forward propagation
% returns: gradients of the cost with respect to each parameter, a struct containing dW1, db1, dW2, db2, etc.

    L = length(parameters);
    m = size(X, 2);
    gradients = cell(1, L);
%     dA = -Y ./ activations{end} + (1 - Y) ./ (1 - activations{end})  ;
%     dZ = dA .* activations{end} .* (1 - activations{end});
    dZ = activations{end} - Y;
    gradients{L}.dW = dZ * activations{L}' / m ;
    gradients{L}.db = sum(dZ, 2) / m ;
    for l = (L-1):-1:1
        dA = parameters{l+1}.W' * dZ;
%         dZ = dA .* (activations{l+1} > 0);
        dZ = dA .* (1 - tanh2(activations{l+1}).^2);
        if l == 1
            A_prev = X;
        else
            A_prev = activations{l};
        end
        gradients{l}.dW = dZ * A_prev' / m ;
        gradients{l}.db = sum(dZ, 2) / m ;
    end
end