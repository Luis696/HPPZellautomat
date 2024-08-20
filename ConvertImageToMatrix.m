% MATLAB code

% Load and convert the image to grayscale
imageMatrix = imread('Rheinahrcampus.jpg');
imageMatrix_gray = rgb2gray(imageMatrix);

% Determine size of the current image
[origRows, origCols] = size(imageMatrix_gray);

% % Determine new size (round up to square root of closest perfect square)
% newSize = ceil(sqrt(origRows*origCols));
newSize = 300;
% If the calculated newSize is less than the length of the row or the column,
% % set it to the larger value between the number of rows and columns
% if newSize < max(origRows, origCols)
%     newSize = max(origRows, origCols);
% end

% Resize image to largest possible square matrix
resizedImage = imresize(imageMatrix_gray, [newSize newSize]);

% Display the image
imshow(resizedImage);

% % Convert the matrix to a double array for saving to a text file
% resizedImage_double = int8(resizedImage);

% Write the matrix to a text file
dlmwrite('resizedImage.txt', resizedImage, 'delimiter', '\t');