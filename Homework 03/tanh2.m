% This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252

% File Name: tanh2.m

function A = tanh2(Z)
% Tanh applies the tanh activation function on the input to get a nonlinear output.
%    Inputs:
%         X: A M x N matrix representing the output of the neurons, which serves as the input of the activation function. M is the number of neurons, and N is the number of examples
%   Outputs:
%         Z: a M x N matrix representing the output after the tanh activation function

    A = 2./(1+exp(-2*Z))-1;
end