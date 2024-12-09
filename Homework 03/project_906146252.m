%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252
% File Name: project_906146252.m

%% Main File
% This is the main script. In this script, you can define network
% architecture and hyperparameters. Running this script will automatically
%load the data from the directory, preprocess it, train your model, test 
% it on the test set and calculate the accuracy of
% your model. And in the end, this training process is visualised through
% plots and matrics, 
clear all; close all; clc;

%% Load Training and Testing Data
[X_train, Y_train, X_test, Y_test] = load_train_and_test_data();

%% Define network architecture and hyperparameters
input_size = size(X_train, 1); % input size of the FNN
output_size = size(Y_train, 1);% output size (number of classes) of the FNN
neurons = 64;  % neurons each layer
numLayer = 5;% No. of hidden layers
lr = 0.01; % learning rate
epochs = 150; % epochs

% to test each case, the variables - no. of hidden layers, learning rate
% and epochs were modified to analyse the performance in each case. 

layer_dims = zeros(1, numLayer + 2);
layer_dims(1) = input_size;
layer_dims(end) = output_size;
for i = 1:numLayer
    layer_dims(i+1) = neurons;
end


%% Train the model
parameters = initialize_parameters(layer_dims);

fprintf('Initialize cost somewhere so that we do not have warnings.\n');

% Train the model using mini-batch gradient descent
m = size(X_train, 2);
batch_size = 64;
num_batches = floor(m / batch_size);
trainLoss = zeros(epochs, 1);
testAccuracy = zeros(epochs, 1);
for i = 1:epochs
    % Shuffle the training data
    indices = randperm(m);
    X_train = X_train(:, indices);
    Y_train = Y_train(:, indices);
    
    % Initialize cost matrix
    cost = zeros(num_batches, batch_size);
    
    % Train on mini-batches
    for j = 1:num_batches
        X_batch = X_train(:, (j-1)*batch_size+1:j*batch_size);
        Y_batch = Y_train(:, (j-1)*batch_size+1:j*batch_size);
        forward_pass = forward_propagation(X_batch, parameters);
        cost(j,:) = compute_cost(forward_pass{end}, Y_batch);
        gradients = backward_propagation(X_batch, Y_batch, parameters, forward_pass);
        parameters = update_parameters(parameters, gradients, lr);
    end
    
    Y_pred = predict(X_test, parameters);
    acc = accuracy(Y_pred, Y_test);
    % Print the cost each epoch
    fprintf('Loss after epoch %d: Training: %f\n', i, norm(cost));
    trainLoss(i) = norm(cost);
    testAccuracy(i) = acc;
end
    

%% Evaluate the model on the test set
fprintf('Test accuracy: %f\n', testAccuracy(end));

%% Visualize the training progress
visualizeHistory(epochs, trainLoss, testAccuracy, lr, numLayer);