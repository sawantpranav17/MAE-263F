%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252

% File Name: softmax.m

function A = softmax(Z)
% softmax applies the softmax function on the input to get the corresponding probability distribution along all possible classes.
%    Inputs:
%         X: A K x N matrix. K is the number of all classes, and N is the number of examples.
%   Outputs:
%         Z: A K x N matrix representing the probability of the distribution. The sum of each row should equal to 1.

    A = exp(Z) ./ sum(exp(Z));
end