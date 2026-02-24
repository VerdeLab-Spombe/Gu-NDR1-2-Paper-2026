% Read the CSV files into tables
win1 = readtable('T-Windows_Ave_Example_1.csv');
pro = readtable('Pro_Ave_example_1.csv');
% 1st win = win(1:57,2:98);
% 2nd win = win(83:122,2:99);
% 3rd win = win(111:164,2:99);
% 4th win = win(166:219,2:99);
% 5th win = win(269:334,2:99);
% 6th win = win(336:401,2:100);

win = win1(1:49,2:100);

pro = pro(1:49,1:99);

data = win;
data = table2array(data);

% Get the size of the data
[rows, cols] = size(data);

% Function to find the closest non-NaN value within a range of 10 rows
function val = find_closest_value(data, row, col, rows)
    max_offset = 10; % Define the maximum offset to search for non-NaN values
    val = NaN; % Initialize the return value as NaN
    for offset = 1:max_offset
        % Check row above
        if row - offset > 0 && ~isnan(data(row - offset, col))
            val = data(row - offset, col);
            return;
        end
        % Check row below
        if row + offset <= rows && ~isnan(data(row + offset, col))
            val = data(row + offset, col);
            return;
        end
    end
    % If no valid value found within 10 rows
    error('More than 10 consecutive NaN values found in column %d at row %d', col, row);
end

% Iterate over each element in the data
for col = 1:cols
    for row = 1:rows
        if isnan(data(row, col))
            data(row, col) = find_closest_value(data, row, col, rows);
        end
    end
end
% Save the modified data back to a CSV file
new_data = array2table(data);
writetable(new_data, '1_layer_Example_1_modified.csv');

win_matrix = table2array(new_data);
pro_matrix = table2array(pro);

% Get the number of rows
[numRows, numCols] = size(win_matrix);

% Set the maximum lag
maxlag = 20;

% Initialize a cell array to store the cross-correlation results
crossCorrelations = cell(numRows, 1);

% Loop through each row and calculate the normalized cross-correlation
for i = 1:numRows
    x = win_matrix(i, :);
    y = pro_matrix(i, :);
    [c, lags] = xcov(x, y, maxlag, 'normalized');
    crossCorrelations{i} = c;
end

% Convert the cell array to a matrix for easier handling
crossCorrelationsMatrix = cell2mat(crossCorrelations);
% Create lag names for the columns
lagNames = compose('Lag_%d', lags);
% Convert the matrix to a table
crossCorrelationsTable = array2table(crossCorrelationsMatrix, 'VariableNames', lagNames);
% Write the table to a CSV file
writetable(crossCorrelationsTable, '1_layer_20_Example_1.csv');
disp('Original data with csv file has been finished.');



% Custom gradient colormap
colors = [
    0 0 0.5;   % Dark blue
    0 0 1;     % Light blue
    0 1 0;     % Green
    1 1 0;     % Yellow
    1 0.5 0;   % Light red
    1 0 0      % Dark red
];
nColors = 256;
customCmap = interp1(linspace(0, 1, size(colors, 1)), colors, linspace(0, 1, nColors));

% Original Heatmap
fig1 = figure;
h1 = heatmap(lags, 1:numRows, crossCorrelationsMatrix);
h1.Title = 'Original Cross-Correlation Heatmap';
h1.XLabel = 'Time Lag (s)';
h1.YLabel = 'Window Number';
h1.ColorLimits = [-0.6, 0.6];
h1.Colormap = customCmap;
h1.GridVisible = 'off';  % Removes the grid lines for a cleaner appearance
h1.CellLabelColor = 'none';  % Hides cell labels to reduce visual clutter

% Save Origina33_layer_20_scr-1103-Cdc42-1-1.fig');
saveas(fig1, '1_layer_20_Example_1.fig');
saveas(fig1, '1_layer_20_Example_1.tif', 'tiff');

% Interpolate data for a smoother appearance

[Xq, Yq] = meshgrid(...
    linspace(1, length(lags), length(lags)*5),...
    linspace(1, numRows, numRows*5));
interpolatedData = interp2(crossCorrelationsMatrix, Xq, Yq, 'cubic');

% Blurred Heatmap
fig2 = figure;
h2 = heatmap(interpolatedData);
%h2.Title = 'Blurred Heatmap';
h2.ColorLimits = [-0.6, 0.6];
h2.Colormap = customCmap;
h2.XLabel = '';  % Remove labels
h2.YLabel = '';
h2.XDisplayLabels = repmat({''}, 1, length(h2.XDisplayLabels));
h2.YDisplayLabels = repmat({''}, 1, length(h2.YDisplayLabels));
h2.GridVisible = 'off';  % Removes the grid lines for a smoother appearance
h2.CellLabelColor = 'none';  % Hides cell labels

% Save Blurred Heatmap1_layer_20_scr-V2-Cdc42-4-2.fig');
saveas(fig2, 'Blurred_1_layer_20_Example_1.fig');
%saveas(fig2, 'Blurred_2_layer_20_scr-1102-Cdc42-1-1.tif', 'tiff');

disp('Both original and blurred heatmaps have been generated with the specified gradient colors.');



