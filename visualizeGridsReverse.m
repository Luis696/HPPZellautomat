clear variables
close all

% Your folder path, replace with your actual folder path
folder_path = 'Grids';

film_speed = 0.01;

% Get a list of all txt files in the folder 
txt_files = dir(fullfile(folder_path, 'grid_*.txt'));

% Create a new figure
figure;

% Define your colormap limits
% These should be defined based on your data
clims = [0 10]; 

% Extract the integer from the filenames and sort
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), {txt_files.name});
[~, sorted_indices] = sort(numbers);

% Loop through each txt file
for k = sorted_indices 
   % Get the file name
   file_name = fullfile(folder_path, txt_files(k).name);

   % Load the data from the file
   matrix_data = load(file_name);

   % Plot data using imagesc
   imagesc(matrix_data, clims); 

   % Add a colorbar for reference
   colorbar;

   % Update the title of the plot with current file name
   title(txt_files(k).name);

   % Pause execution for a short period of time to allow plot to update
   pause(film_speed);
end

% Your folder path, replace with your actual folder path
folder_path = 'Grids';

% Create a new figure
figure;

first_encryption = fullfile(folder_path, txt_files(1).name);
first_encryption = load(first_encryption);

last_decryption = fullfile(folder_path, txt_files(end).name);
last_decryption = load(last_decryption);

figure()
title("Decription - Encryption")
imagesc(first_encryption-last_decryption)
colorbar
