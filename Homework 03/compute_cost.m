%This code and the helper functions/files have been referenced
%from the class resources on Bruinlearn

%Name: Pranav Sawant
%UID: 906146252
%%File Name: compute_cost.m

function cost = compute_cost(AL, Y)
%The given function calculates the cross-entropy loss, also knows as the
%log loss, given the predicted values AL and the true values Y
%Input: AL: Final predicted values, Y:Ground Truth Labels
%Output: cost: entropy loss
    cost = -sum(Y .* log(AL));
end