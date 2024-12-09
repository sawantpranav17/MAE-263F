% This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252
% File Name: visualizeHistory.m

function visualizeHistory(epochs, trainLoss, testAccuracy, lr, numLayer)
% visualizeHistory: a function that plots and saves the training loss and testing accuracy versus time under specific hyperparameters.
%   Inputs:
%       epochs: a integer shows the number of training epochs
%       trainLoss: an epochs*1 vector, where ith element in the vector corresponds to the training loss after ith epochs training
%       testAccuracy: an epochs*1 vector, where ith element in the vector corresponds to the testing accuracy after ith epochs training
%       lr: the value of learning rate
%       numLayer: the value of the number of hidden layers
%   Outputs:
%       This function has no outputs

h1 = figure(1);
subplot(1, 2, 1)
plot(1:epochs, trainLoss)
xlabel('Epochs');
ylabel('Training Loss');
subplot(1, 2, 2)
plot(1:epochs, testAccuracy)
xlabel('Epochs');
ylabel('Testing Accuracy');


sgtitle(sprintf("Epochs: %d; learning rate: %g; number of hidden layers: %d", ...
    epochs, lr, numLayer));

filename = sprintf("model_%g_%d_%d.png", lr, numLayer, epochs);
saveas(h1, filename);


end
