%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252
% File Name: load_train_and_test_data.m

function [X_train, Y_train, X_test, Y_test] = load_train_and_test_data()
% Reads training and testing data. The images are flattened out and normalized, labels are converted to one-hot encoded vector.
%     Inputs:
%         (None)
%     Outputs:
%     X_train: a  784 x N matrix representing the training images, where each column corresponds to a single image. N is the number of examples(images), which is 60000. The pixel intensities are normalized, holding values from 0 to 1.
%     Y_train: a K x N matrix representing the labels of the training images, where each column corresponds to a one-hot encoding of a single image. K is the number of classes, which is 10. N is the number of examples, 60000.
%      X_test: a 784 x N_t matrix representing the testing images, where each column corresponds to a single image. The pixel intensities are normalized, holding values from 0 to 1. N_t is the number of test examples, which is 100000.
%      Y_test: a K x N_t matrix representing the labels of the testing images, where each column corresponds to a one-hot encoding of a single image.

    % Load training data
    load('train_images.mat');
    load('train_labels.mat');
    pixel = reshape(pixel, [size(pixel, 1)*size(pixel, 2), size(pixel, 3)]) / 255;
    X_train = pixel;
    % Convert y to one-hot encoding
    Y_train = zeros(10, length(label));
    for i = 1:length(label)
        Y_train(label(i)+1, i) = 1;
    end
    
    % Load testing data
    load("test_images.mat");
    load("test_labels.mat");
    pixel = reshape(pixel, [size(pixel, 1)*size(pixel, 2), size(pixel, 3)]) / 255;
    X_test = pixel;
    Y_test = zeros(10, length(label));
    for i = 1:length(label)
        Y_test(label(i)+1, i) = 1;
    end
    
end