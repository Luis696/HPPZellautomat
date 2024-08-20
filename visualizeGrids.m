clear variables
close all

% Your folder path, replace with your actual folder path
folder_path = 'Grids';

film_speed = 0.05;

% Get a list of all txt files in the folder 
txt_files = dir(fullfile(folder_path, 'encrypting_grid_*.txt'));

% Create a new figure
figure;

% Define your colormap limits
% These should be defined based on your data
clims = [0 16]; 

% Extract the integer from the filenames and sort
numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), {txt_files.name});
[~, sorted_indices] = sort(numbers);
% Create a new figure
    figure;
% Loop through each txt file
for k = sorted_indices 
    % figure;
    % Get the file name
   file_name = fullfile(folder_path, txt_files(k).name);

   % Load the data from the file
   matrix_data = load(file_name);
    
   % Plot data using imagesc
   imagesc(matrix_data, clims); 
   
     % Erhalte die Größe des Arrays
    [numRows, numCols] = size(matrix_data);

    % Schreibe die Zahlen auf das Bild
    for row = 1:numRows
        for col = 1:numCols
            num = matrix_data(row, col);
            text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
        end
    end

   % Add a colorbar for reference
   colorbar;

   % Update the title of the plot with current file name
   title(txt_files(k).name);

   % Pause execution for a short period of time to allow plot to update
   % pause(film_speed);
   waitforbuttonpress;
end


first_encryption = fullfile(folder_path, txt_files(1).name);
first_encryption = load(first_encryption);

% Your folder path, replace with your actual folde
% r path
folder_path = 'Grids';

% Get a list of all txt files in the folder 
txt_files = dir(fullfile(folder_path, 'decrypting_grid_*.txt'));

% Create a new figure
figure;



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
   

     % Erhalte die Größe des Arrays
    [numRows, numCols] = size(matrix_data);

    % Schreibe die Zahlen auf das Bild
    for row = 1:numRows
        for col = 1:numCols
            num = matrix_data(row, col);
            text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
        end
    end
    

   % Add a colorbar for reference
   colorbar;

   % Update the title of the plot with current file name
   title(txt_files(k).name);


   % Add a colorbar for reference
   colorbar;

   % Update the title of the plot with current file name
   title(txt_files(k).name);

   % Pause execution for a short period of time to allow plot to update
   waitforbuttonpress;
end

last_decription = load(file_name);



% check for failures in particle movement:

figure()
Decrypt_Map= [0, 4, 8, 12, 1, 10, 9, 13, 2, 6, 5, 14, 3, 7, 11, 15];
% Erhalte die Größe des Arrays
[numRows, numCols] = size(last_decription);
diff_image = zeros(numRows,numCols);
% Schreibe die Zahlen auf das Bild
for row = 1:numRows
    for col = 1:numCols
        diff = first_encryption(row,col) - Decrypt_Map(last_decription(row, col)+1);
        diff_image(row,col) = diff;
    end
end

imagesc(diff_image)
title("difference between encryption and decryption")


figure()
subplot(2,1,1)
imagesc(first_encryption)
subtitle("encryption") 
% Erhalte die Größe des Arrays
[numRows, numCols] = size(first_encryption);
 
% Schreibe die Zahlen auf das Bild
for row = 1:numRows
    for col = 1:numCols
        num = first_encryption(row, col);
        if diff_image(row, col) ~= 0
            text(col, row, num2str(num), 'Color', 'red', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        else
            % text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
        end
    end
end


subplot(2,1,2)
imagesc(last_decription)
subtitle("decryption")
% Erhalte die Größe des Arrays
[numRows, numCols] = size(last_decription);

% Schreibe die Zahlen auf das Bild
for row = 1:numRows
    for col = 1:numCols
        num = last_decription(row, col);
        if diff_image(row, col) ~= 0
            text(col, row, num2str(num), 'Color', 'red', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        else
            % text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
        end
    end
end


