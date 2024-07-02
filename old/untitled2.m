% Find peaks to identify each period of the signal
data_norm = (AP_Amp)./max(AP_Amp);
data_norm_smooth = smoothdata(data_norm, 'sgolay', 5);
[pks,locs,w,p] = findpeaks(data_norm_smooth, 'MinPeakDistance',300, 'MinPeakProminence', 0.01);

% Visualize the peaks
h = figure; 
h.Position = [224, 537, 1215, 565]; 
plot(data_norm_smooth, 'Color', [0.8, 0.8, 0.8]); % Plot all signals in gray
hold on; grid on;
plot(locs, pks, 'r*'); % Plot peaks

% Iterate over each period to stack them and find the max length
max_length = 0;
for iter1 = 1:length(locs)-1
    temp_data = data_norm_smooth(locs(iter1):locs(iter1+1));
    length_temp = length(temp_data);
    if length_temp > max_length
        max_length = length_temp;
    end
end

% Preallocate a matrix to store aligned signals
gatingSignalSet = nan(length(locs)-1, max_length);

% Align each signal period and store it in the matrix
for iter1 = 1:length(locs)-1
    temp_data = data_norm_smooth(locs(iter1):locs(iter1+1));
    gatingSignalSet(iter1, 1:length(temp_data)) = temp_data;
    % Plot each period signal in gray
    plot([0:length(temp_data)-1].*x_step_second, temp_data, 'Color', [0.5, 0.5, 0.5]);
end

% Calculate the average signal from the aligned signals
gatingSignalSetAvg = nanmean(gatingSignalSet, 1);
gatingSignalSetAvg_min = min(gatingSignalSetAvg);
gatingSignalSetAvg_max = max(gatingSignalSetAvg);

% Normalize the average signal
gatingSignalSetAvg_normalized = (gatingSignalSetAvg - gatingSignalSetAvg_min) / (gatingSignalSetAvg_max - gatingSignalSetAvg_min);

% Calculate the time step in seconds
x_step_second = mean(diff(AP_time));

% Generate the time data for plotting
XData_second = [0:size(gatingSignalSetAvg_normalized, 2)-1].*x_step_second;

% Define the phase area to highlight
startPhase = 25; % Example starting phase percentage
endPhase = 75;   % Example ending phase percentage

% Convert phase percentages to indices
startIndex = round((startPhase/100) * length(gatingSignalSetAvg_normalized));
endIndex = round((endPhase/100) * length(gatingSignalSetAvg_normalized));

% Plot the normalized average signal in black dashed line
plot(XData_second, gatingSignalSetAvg_normalized, 'k--'), 
xlabel("Time [s]"), ylabel("Normalized Amplitude"), xlim([0, max(XData_second)]);

% Add yellow patch for the specified phase range
ylimits = get(gca, 'YLim');
patch([XData_second(startIndex) XData_second(endIndex) XData_second(endIndex) XData_second(startIndex)], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], 'yellow', 'FaceAlpha', 0.3);

% Ensure the average signal line is on top
uistack(findall(gca, 'Type', 'line'), 'top');

legend('Periodic signals', 'Peaks', 'Average Signal');

% 축의 위치를 가져옵니다.
axPos = get(gca, 'Position'); % [left bottom width height]

% 현재 축의 x 한계와 y 한계를 가져옵니다.
xlimits = xlim(gca);
ylimits = ylim(gca);

% 데이터 좌표에서 정규화된 그림 좌표로 변환합니다.
xStartNorm = (XData_second(startIndex) - xlimits(1)) / (xlimits(2) - xlimits(1));
xEndNorm = (XData_second(endIndex) - xlimits(1)) / (xlimits(2) - xlimits(1));
yPosNorm = (0.5 - ylimits(1)) / (ylimits(2) - ylimits(1));

% 축 위치를 기반으로 최종 정규화된 그림 좌표를 계산합니다.
xStartNorm = axPos(1) + xStartNorm * axPos(3);
xEndNorm = axPos(1) + xEndNorm * axPos(3);
yPosNorm = axPos(2) + yPosNorm * axPos(4);

% 이제 정규화된 좌표를 사용하여 double arrow를 추가합니다.
annotation('doublearrow', [xStartNorm xEndNorm], [yPosNorm yPosNorm], 'Color', 'k');

% Add yellow patch for the specified phase range
ylimits = get(gca, 'YLim');
patch([XData_second(startIndex) XData_second(endIndex) XData_second(endIndex) XData_second(startIndex)], [ylimits(1) ylimits(1) ylimits(2) ylimits(2)], 'yellow', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Add double-headed arrow
yPos = 0.5 * (ylimits(1) + ylimits(2)); % Position for the arrow in data units
annotation('doublearrow', [XData_second(startIndex) XData_second(endIndex)], [yPos yPos], 'Color', 'k', 'Head1Length', 10, 'Head1Width', 10, 'Head2Length', 10, 'Head2Width', 10);

% Ensure the average signal line and peaks are on top
uistack(findall(gca, 'Type', 'line'), 'top');

% Optional: Add text label for the double arrow if needed
text(mean([XData_second(startIndex) XData_second(endIndex)]), yPos, sprintf('Time: %.2f s', XData_second(endIndex)-XData_second(startIndex)), 'HorizontalAlignment', 'center', 'BackgroundColor', 'none', 'EdgeColor', 'none');

% Redraw the plot to ensure the lines are above the shaded area
hold off;

% Apply legend if needed
legend('Periodic signals', 'Peaks', 'Average Signal', 'Time interval');

