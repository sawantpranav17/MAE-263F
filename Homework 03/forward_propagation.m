%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252
% File Name: forward_propagation.m

function activations = forward_propagation(X, parameters)
% Computes the output of a feedforward neural network given input data and learned parameters
% X: input data, shape (input size, number of examples)
% parameters: learned parameters, a struct containing W1, b1, W2, b2, etc.
% returns: array of activations at each layer, including input and output layers

    L = length(parameters);
    A = X;
    activations = cell(1, L+1);
    activations{1} = X;
    for l = 1:L
        Z = parameters{l}.W*A + parameters{l}.b;
%         Z = linear_forward(A, parameters{l}.W, parameters{l}.b);
        if l == L
            A = softmax(Z);
        else
%             A = sigmoid(Z);
            A = tanh2(Z);
        end
        activations{l+1} = A;
    end
end