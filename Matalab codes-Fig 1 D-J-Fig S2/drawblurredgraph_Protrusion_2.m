
% Read the CSV files into tables
data1 = readtable('Pro_Ave_example_1.csv');
data = data1(:, 1:91);
%data = data1;
crossCorrelationsMatrix = table2array(data);
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
%Sequence = 0:10:980; 
xLabels = 0:10:900;
% Original x values from 0 to 900 with an interval of 10
customXLabels = arrayfun(@(x) num2str(x), 0:100:900, 'UniformOutput', false);

[numRows, numCols] = size(crossCorrelationsMatrix);

% Define the custom y-axis labels
maxYValue = ceil(numRows / 10) * 10; % Round up to the next multiple of 10
yLabels = 0:10:maxYValue; % Custom y values from 0 to maxYValue with an interval of 10
customYLabels = arrayfun(@(y) num2str(y), yLabels, 'UniformOutput', false);



% Automatically set color limits to the data range
dataMin = min(crossCorrelationsMatrix(:));
dataMax = max(crossCorrelationsMatrix(:));
h.ColorLimits = [dataMin, dataMax];



% Original Heatmap
fig1 = figure;
h1 = heatmap(xLabels, 1:numRows, crossCorrelationsMatrix);
h1.Title = 'Original Heatmap';
h1.XLabel = 'Time (s)';
h1.YLabel = 'Window Number';
h1.ColorLimits = [-10, 10];
%h1.ColorLimits = [dataMin, dataMax];
h1.Colormap = customCmap;
h1.GridVisible = 'off';  % Removes the grid lines for a cleaner appearance
h1.CellLabelColor = 'none';  % Hides cell labels to reduce visual clutter
h1.XData = xLabels;
h1.XDisplayLabels = repmat({''}, size(xLabels)); % Empty labels for all positions
h1.XDisplayLabels(1:10:end) = customXLabels; % Set custom labels at every 200 interval


% Set the YData and YDisplayLabels properties
h1.YData = 1:numRows; % Assuming these are the original row indices
h1.YDisplayLabels = repmat({''}, size(1:numRows)); % Empty labels for all positions
% Set custom labels at every 10 interval within the available range
for i = 1:length(yLabels)
    index = yLabels(i) + 1; % Convert yLabel value to index (MATLAB is 1-based)
    if index <= numRows
        h1.YDisplayLabels{index} = customYLabels{i};
    end
end



% Save Original Heatmap
saveas(fig1, 'Adjusted_Pro_Ave_example_1.fig');
%saveas(fig1, '3_layer_Windows_1102-cdc42-1-1-1.tif', 'tiff');

% Interpolate data for a smoother appearance

[Xq, Yq] = meshgrid(...
    linspace(1, length(xLabels), length(xLabels)*5),...
    linspace(1, numRows, numRows*5));
interpolatedData = interp2(crossCorrelationsMatrix, Xq, Yq, 'cubic');

% Blurred Heatmap
fig2 = figure;
h2 = heatmap(interpolatedData);
%h2.Title = 'Blurred Heatmap';
h2.ColorLimits = [-10, 10];
%h2.ColorLimits = [dataMin, dataMax];
h2.Colormap = customCmap;
h2.XLabel = '';  % Remove labels
h2.YLabel = '';
h2.XDisplayLabels = repmat({''}, 1, length(h2.XDisplayLabels));
h2.YDisplayLabels = repmat({''}, 1, length(h2.YDisplayLabels));
h2.GridVisible = 'off';  % Removes the grid lines for a smoother appearance
h2.CellLabelColor = 'none';  % Hides cell labels


% Save Blurred Heatmap
saveas(fig2, 'Adjusted_Blurred_Pro_Ave_example_1.fig');
%saveas(fig2, 'Adjusted_Blurred_1_layer_20_1103-RhoA-2-7.tif', 'tiff');

disp('Both original and blurred heatmaps have been generated with the specified gradient colors.');