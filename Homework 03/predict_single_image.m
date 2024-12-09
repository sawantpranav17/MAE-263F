%This code and the helper functions/files have been referenced from the class resources on Bruinlearn

% Name: Pranav Sawant
% UID: 906146252

% File Name: predict_single_image.m

function predict_single_image(parameters)
% Visualizes the input image and plots a bar graph of probabilities for
% each class
% Inputs:
%         parameters: a struct containing the weights and biases (W1, b1,
%         W2, b2, etc.)
% Output: 
%         (None)

    % Load data
    load('test_images.mat');
    index = randi(length(pixel));
    X = reshape( pixel(:,:,index), [size(pixel, 1)*size(pixel, 2), 1]);

    forward = forward_propagation(X, parameters);
    probabilities = forward{end};

    figure();
    subplot(1, 2, 1)
    imshow(pixel(:,:,index))
    title('Input Image')
    subplot(1, 2, 2)
    bar(0:9, probabilities)
    xlabel('Classes');
    ylabel('Probability');
    title('Probability Distribution')

end