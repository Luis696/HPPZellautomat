clear variables
close all
time = 0.0001;
% Your folder path, replace with your actual folder path
folder_path = 'Grids';


% Get a list of all txt files in the folder 
txt_files = dir(fullfile(folder_path, 'encrypting_message_*.txt'));

% Create a new figure
figure();

% Define your colormap limits
% These should be defined based on your data
 

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
    % imshow(uint8(matrix_data)); 
    imagesc(matrix_data);

    %  % Erhalte die Größe des Arrays
    % [numRows, numCols] = size(matrix_data);
    % 
    % % Schreibe die Zahlen auf das Bild
    % for row = 1:numRows
    %     for col = 1:numCols
    %         num = matrix_data(row, col);
    %         text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
    %     end
    % end

   % Add a colorbar for reference
   % colorbar;

   % Update the title of the plot with current file name
   title(txt_files(k).name);

   % Pause execution for a short period of time to allow plot to update
   % pause(film_speed);
   pause(time);
end

first_encryption = fullfile(folder_path, txt_files(1).name);
first_encryption = load(first_encryption);

% Your folder path, replace with your actual folde
% r path
folder_path = 'Grids';

% Get a list of all txt files in the folder 
txt_files = dir(fullfile(folder_path, 'decrypting_message_*.txt'));

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
   imagesc(matrix_data); 
   % imshow(uint8(matrix_data));
    % 
    %  % Erhalte die Größe des Arrays
    % [numRows, numCols] = size(matrix_data);

    % % Schreibe die Zahlen auf das Bild
    % for row = 1:numRows
    %     for col = 1:numCols
    %         num = matrix_data(row, col);
    %         text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
    %     end
    % end

   % % Add a colorbar for reference
   % colorbar;

   % Update the title of the plot with current file name
   title(txt_files(k).name);


   

   % Update the title of the plot with current file name
   title(txt_files(k).name);

   % Pause execution for a short period of time to allow plot to update
   pause(time);
end

last_decription = load(file_name);


figure()
subplot(2,1,1)
subtitle("encryption")
imagesc(first_encryption)
% Erhalte die Größe des Arrays
[numRows, numCols] = size(first_encryption);

% Schreibe die Zahlen auf das Bild
for row = 1:numRows
    for col = 1:numCols
        num = first_encryption(row, col);
        text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
    end
end


subplot(2,1,2)
subtitle("decryption")
imagesc(last_decription)
% Erhalte die Größe des Arrays
[numRows, numCols] = size(last_decription);

% Schreibe die Zahlen auf das Bild
for row = 1:numRows
    for col = 1:numCols
        num = last_decription(row, col);
        text(col, row, num2str(num), 'Color', 'black', 'HorizontalAlignment', 'center');
    end
end
