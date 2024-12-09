%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252
% File Name: linear_forward.m

function Z = linear_forward(A_prev, W, b)
%Helper function for implementing forward path. For each hidden layer A is
%calculated and fed into the next layer as A_prev
%Inputs: A_prev: A caculated from the previous layer
%b: Bias term, W: weigth of the connection
    Z = W * A_prev + b;
end