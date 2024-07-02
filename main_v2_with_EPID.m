clc
clear
close all
fclose all
%%
% ???? ??? ????
% filename = "./QA_2403_25t75.txt";
% filename = "./240326_GatingQA_Signal_OfflineReview.txt"
filename = "\\10.10.37.115\Research\0.wjcheon\00_Research\07_GatingQA @ynkang\EPID video\trial_3\signal_240624_3.txt"
filename = "QA_2403_25t75.txt"
st_phase = 25;
end_phase = 75;
filenameFrameSum = "\\10.10.37.115\Research\0.wjcheon\00_Research\07_GatingQA @ynkang\EPID video\trial_3\frame_sums_trail3.csv"
filenameFrameSum = "frame_sums2.csv"

% ???? ??? ????
leftRightName = "Left-Right";
anteriorPosteriorName = "Anterior-Posterior";
headFeetName = "Head-Feet";
phaseChannelName = "Phase Channel";
beamFlagName = "Beam Enable/Disable Moments";
%% Read the log files  
% ???? ????
fileID = fopen(filename, 'r');
if fileID == -1
    error('File cannot be opened: %s', filename);
end

try
    while ~feof(fileID)
        lineTemp = strtrim(fgetl(fileID)); % Trim leading/trailing whitespaces
        if strcmp(lineTemp, anteriorPosteriorName)
            % Skip 3 lines
            for i = 1:3
                fgetl(fileID);
            end
            % ap
            AP_time = [];
            AP_Amp = [];
            while true
                lineTemp = strtrim(fgetl(fileID));
                if isempty(lineTemp) || strcmp(lineTemp, "--------------------------------------------------------------------------------------")
                    break; % Exit loop if empty line or next section starts
                end
                [C, ~] = strsplit(lineTemp, "\t");
                AP_time(end+1) = str2double(C{1});
                AP_Amp(end+1) = str2double(C{2});
            end
            dataStruct.LR = [AP_time; AP_Amp];
            % You can add more sections here following the same pattern
        end
    end
catch
    warning('Error while reading the file:');
end

% ???? ???
fclose(fileID);
% h0 = figure, 


%

% ???? ????
fileID = fopen(filename, 'r');
if fileID == -1
    error('File cannot be opened: %s', filename);
end

try
    while ~feof(fileID)
        lineTemp = strtrim(fgetl(fileID)); % Trim leading/trailing whitespaces
        if strcmp(lineTemp, phaseChannelName)
            % Skip 3 lines
            for i = 1:3
                fgetl(fileID);
            end
            % Phase
            phase_time = [];
            phase_Amp = [];
            while true
                lineTemp = strtrim(fgetl(fileID));
                if isempty(lineTemp) || strcmp(lineTemp, "--------------------------------------------------------------------------------------")
                    break; % Exit loop if empty line or next section starts
                end
                [C, ~] = strsplit(lineTemp, "\t");
                phase_time(end+1) = str2double(C{1});
                phase_Amp(end+1) = str2double(C{2});
            end
            dataStruct.LR = [phase_time; phase_Amp];
            % You can add more sections here following the same pattern
        end
    end
catch
    warning('Error while reading the file:');
end

% ???? ???
fclose(fileID);
% figure, plot(phase_time, phase_Amp)


%
% ???? ????
fileID = fopen(filename, 'r');
if fileID == -1
    error('File cannot be opened: %s', filename);
end

try
    while ~feof(fileID)
        lineTemp = strtrim(fgetl(fileID)); % Trim leading/trailing whitespaces
        if strcmp(lineTemp, beamFlagName)
            % Skip 3 lines
            for i = 1:3
                fgetl(fileID);
            end
            % Read Left-Right data
            beam_time = [];
            beam_Amp = [];
            while true
                lineTemp = strtrim(fgetl(fileID));
                if isempty(lineTemp) || strcmp(lineTemp, "--------------------------------------------------------------------------------------")
                    break; % Exit loop if empty line or next section starts
                end
                [C, ~] = strsplit(lineTemp, "\t");
                beam_time(end+1) = str2double(C{1});
                beam_Amp(end+1) = str2double(C{2});
            end
            dataStruct.beamTrigger = [beam_time; beam_Amp];
            % You can add more sections here following the same pattern
        end
    end
catch
    warning('Error while reading the file:');
end

% ???? ???
fclose(fileID);
% figure, plot(beam_time, beam_Amp, 'r*')

%
% Initialize beam status as off for all AP_time points
beam_status = zeros(size(AP_time)); 
% Loop through beam_Amp to find transitions, ensuring we don't go out of bounds
for i = 1:(length(beam_Amp)-1)
    if beam_Amp(i) == 1 % If the beam is on at this point
        % Attempt to find the next off point, safely handling the end of array
        offIndex = find(beam_Amp((i+1):end) == 0, 1, 'first');
        if isempty(offIndex)
            % If the beam doesn't turn off again, consider it on until the end of beam_time
            offTime = beam_time(end);
        else
            % Correctly calculate offTime without exceeding array bounds
            offTime = beam_time(i + offIndex);
        end
        
        % Find the range in AP_time that corresponds to the current on period
        indices = find(AP_time >= beam_time(i) & AP_time <= offTime);
        
        % Check if indices are not empty to avoid setting an empty slice
        if ~isempty(indices)
            % Set beam status to 1 (on) for these indices
            beam_status(indices) = 1;
        end
    end
end


%%
% Video read 
%
close all
FrameSumData = xlsread(filenameFrameSum );
times = FrameSumData(:,1);
frameSums = FrameSumData(:,2);

%Pre-processing
% minVal = min(frameSums);
% maxVal = max(frameSums);
% 
% % �ּҰ��� �ִ밪�� ������
% frameSums = maxVal + minVal - frameSums;
% 
% % ����ȭ: �ִ밪�� 1�� �ǵ��� ����
% frameSums = frameSums / max(frameSums);


frameSums_Binary = frameSums>5.8 *10^6;
% frameSums_Binary = frameSums>0.8;
% figure, plot(times, frameSums)
% figure, plot(times, frameSums_Binary)

frameSums_Bianry_diff = diff(frameSums_Binary);

% figure, plot(times(1:end-1), frameSums_Bianry_diff )

frameSums_Bianry_diff_1 = find(frameSums_Bianry_diff==1);
frameSums_Bianry_diff_0 = find(frameSums_Bianry_diff==-1);
frameSums_Bianry_diff_0(1) = [];
frameSums_Bianry_diff_1(end) = [];
% frameSums_Bianry_diff_1(end) = [];

index_diff = frameSums_Bianry_diff_0-frameSums_Bianry_diff_1;

frameRate = mean(diff(times));
mean(index_diff.*frameRate)
std(index_diff.*frameRate)




% Now, beam_status contains 1s at indices where the beam is on and 0s where it is off
%
figure, hold on
% p1 = plot(AP_time, beam_status);
p2 = plot(AP_time, AP_Amp), title('AP'), hold on
p3 = plot(AP_time, beam_status.*0.07), title('Beam on/off')
p4 = plot(times, frameSums_Binary*0.097)
legend([p2, p3, p4], {'amp1', 'BEAM TRIGGER', 'epdi'})
grid on


time_pad = [0:99]*0.04
time_new = [time_pad'; times]

frameSums_pad = zeros(100,1);
frameSums_new = [frameSums_pad; frameSums]

interp_frameSums = interp1(time_new, frameSums_new, AP_time)


%% Generates averaged periodic signal
% AP_Amp = AP_Amp(1:3740)
% AP_time = AP_time(1:3740)

% beam_status(3726:end) = 0;



% find peaks 
AP_Amp_min = min(AP_Amp(:));
AP_Amp_max = max(AP_Amp(:));
AP_Amp_norm = (AP_Amp-AP_Amp_min) / (AP_Amp_max - AP_Amp_min);

% data_norm = (AP_Amp)./max(AP_Amp);
data_norm = AP_Amp_norm;
data_norm_smooth = smoothdata(data_norm, 5);
[pks,locs,w,p] = findpeaks(data_norm, 'MinPeakDistance',300, 'MinPeakProminence', 0.01)

% visulaization the peaks
h=figure, 
h.Position = [223.7924528301887,537.1509433962264,1215.207547169811,564.8490566037736] 
plot(data_norm), hold on, grid on 
plot(locs, pks, 'r*')

max_length = 0
for iter1=1:length(locs)-1
    temp_data = data_norm_smooth(locs(iter1):locs(iter1+1));
    length_temp = length(temp_data);
    if length_temp > max_length
        max_length = length_temp;
    end
end

gatingSignalSet = nan(length(locs)-1, max_length);
epidSignalSet = nan(length(locs)-1, max_length);
for iter1=1:length(locs)-1
    temp_data = data_norm_smooth(locs(iter1):locs(iter1+1));
    temp_data_epid = interp_frameSums(locs(iter1):locs(iter1+1));
    gatingSignalSet(iter1, 1:length(temp_data)) = temp_data;
    epidSignalSet(iter1, 1:length(temp_data)) = temp_data_epid;
end


gatingSignalSet_min = min(gatingSignalSet(:));
gatingSignalSet_max = max(gatingSignalSet(:));
gatingSignalSet_normalized = (gatingSignalSet - gatingSignalSet_min) / (gatingSignalSet_max - gatingSignalSet_min);

gatingSignalSetAvg = nanmean(gatingSignalSet_normalized, 1)
% gatingSignalSetAvg_min = min(gatingSignalSetAvg);
% gatingSignalSetAvg_max = max(gatingSignalSetAvg);
% gatingSignalSetAvg_normalized = (gatingSignalSetAvg - gatingSignalSetAvg_min) / (gatingSignalSetAvg_max - gatingSignalSetAvg_min);
gatingSignalSetAvg_normalized = gatingSignalSetAvg;


x_step_second = mean(diff(AP_time))
XData_second = [0:size(gatingSignalSetAvg_normalized, 2)-1].*x_step_second;
% XData_second = AP_time(1:3726)
% 
figure, plot(XData_second, gatingSignalSetAvg_normalized, 'k--'), grid on, xlabel("Time [s]"), ylabel("Normalized amplitude"), xlim([0, max(XData_second)]), 
hold on

epidSignalSet_avg = mean(epidSignalSet(2:end-1,1:end), 1)
figure, plot(XData_second, epidSignalSet_avg, 'k--'), grid on, xlabel("Time [s]"), ylabel("AVG. EPID signal")

%%

% Close any open figures
close all;

% Assuming the gatingSignalSetAvg_normalized and signalInterval are already defined
originalIndicesPhase = linspace(0, 100, length(gatingSignalSetAvg_normalized));

% Calculate interpolated values
interpolatedValues_25 = interp1(originalIndicesPhase, gatingSignalSetAvg_normalized, st_phase, 'linear');
interpolatedValues_75 = interp1(originalIndicesPhase, gatingSignalSetAvg_normalized, end_phase, 'linear');

% Assuming signalInterval is defined
signalInterval = x_step_second; % Replace with the actual interval value

% Calculate the times at the 25th and 75th percentiles
timeAt25 = signalInterval * (st_phase/100 * length(gatingSignalSetAvg_normalized) - 1);
timeAt75 = signalInterval * (end_phase/100 * length(gatingSignalSetAvg_normalized) - 1);
timeDifference = timeAt75 - timeAt25;
% Close any open figures
close all;

% Assuming the gatingSignalSetAvg_normalized and signalInterval are already defined
originalIndicesPhase = linspace(0, 100, length(gatingSignalSetAvg_normalized));

% Calculate interpolated values
interpolatedValues_25 = interp1(originalIndicesPhase, gatingSignalSetAvg_normalized, st_phase, 'linear');
interpolatedValues_75 = interp1(originalIndicesPhase, gatingSignalSetAvg_normalized, end_phase, 'linear');

% Assuming signalInterval is defined
signalInterval = x_step_second; % Replace with the actual interval value

% Calculate the times at the 25th and 75th percentiles
timeAt25 = signalInterval * (st_phase/100 * length(gatingSignalSetAvg_normalized) - 1);
timeAt75 = signalInterval * (end_phase/100 * length(gatingSignalSetAvg_normalized) - 1);
timeDifference = timeAt75 - timeAt25;

% Create a new figure and plot
hhh0 =figure;
h1 = plot(originalIndicesPhase, gatingSignalSetAvg_normalized, 'k--', "LineWidth", 2);
hold on; % Keep the plot active for further plotting
grid on;
xlabel("Phase (%)");
ylabel("Normalized Amplitude");
title('Average of periodic signal');

% Define the y-limits
yl = ylim;

% Create the yellow patch
patch([st_phase end_phase end_phase st_phase], [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Get the current axis limits
xLimits = xlim;
yLimits = ylim;

axPos = get(gca, 'Position');

% Convert the X data points of the phase to normalized figure units
x_fig_start = axPos(1) + (st_phase/100 * axPos(3));
x_fig_end = axPos(1) + (end_phase/100 * axPos(3));

% The Y position for the annotation arrow in normalized figure units
% We use the middle of the Y-axis range for the vertical position
y_fig = axPos(2) + (0.5 - yLimits(1))/diff(yLimits) * axPos(4);

% Add a double-headed arrow annotation at y = 0.5 within the yellow shaded region
annotation('doublearrow', [x_fig_start, x_fig_end], [y_fig, y_fig], 'Color', 'black');

% Assuming beamOnTime is defined
% Assuming beamOnTime is defined
text(50, 0.6, sprintf("Time between %d%% to %d%% phases: %1.2f s", st_phase, end_phase, timeDifference), 'HorizontalAlignment', 'center', 'BackgroundColor', 'white');

% Remove the hold on state
% hold off;
plot(originalIndicesPhase, gatingSignalSet_normalized', "Color", [0.85, 0.85, 0.85], "LineStyle", "--")
hold off

% Reapply the legend
legend({sprintf("Time between %d%% to %d%% phases", st_phase, end_phase), "Average of periodic signal"});
% saveas(hhh0,  "????? ???? ?????????_1 (100ms ????).png")
%%

% close all 
% 
% AP_time_new = AP_time(1:3740);
% beam_status = beam_status(1:3740);

% Find changes in the beam status to identify beam on (1) and off (0) edges
statusChanges = diff([0, beam_status, 0]); % Pad with 0 to capture change at start and end
onEdges = find(statusChanges == 1);
offEdges = find(statusChanges == -1);
% offEdges(end) = 2068;

% Initialize array to store durations of 'beam on' states
beamOnDurations = zeros(length(onEdges), 1);

% Calculate the duration of each 'beam on' state
for i = 1:length(onEdges)
    % Duration is the difference in AP_time between off edge and on edge
    beamOnDurations(i) = AP_time(offEdges(i)) - AP_time(onEdges(i));
end

% Display the results
for i = 1:length(beamOnDurations)
    fprintf('Beam On Group %d: Duration = %f seconds\n', i, beamOnDurations(i));
end

% Remove the maximum and minimum values
durationsWithoutExtremes = beamOnDurations(beamOnDurations ~= max(beamOnDurations) & beamOnDurations ~= min(beamOnDurations));

% Calculate the average of the remaining durations
if isempty(durationsWithoutExtremes)
    avgDuration = 0; % This means there were only two groups or all durations were the same
    disp('Not enough data to exclude max and min values.');
else
    avgDuration = mean(durationsWithoutExtremes);
end

% Display the average duration
fprintf('Average Beam On Duration (excluding extremes): %f seconds\n', avgDuration);

% Create the plot
hh = figure; 
ax = axes('Parent', hh);
plot(ax, AP_time, beam_status, 'k--');
xlabel('Time (s)');
ylabel('Beam status');
ylim([0 1.1]); % Adjust the y-axis limits to give some space for annotations

% Y?? ???? 1?? ?????? ????? ????????? ????
hold on;
for i = 1:length(onEdges)
    patch([AP_time(onEdges(i)) AP_time(offEdges(i)) AP_time(offEdges(i)) AP_time(onEdges(i))], [0 0 1 1], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
end

% ?? ????? ?????????? doublearrow?? ????? ?? ?????? ???
for i = 1:length(onEdges)
    % Get the position of the plot in normalized units
    ax_pos = get(ax, 'Position');
    
    % Convert start and end times of the beam on to normalized units
    norm_start = (AP_time(onEdges(i)) - ax.XLim(1)) / range(ax.XLim);
    norm_end = (AP_time(offEdges(i)) - ax.XLim(1)) / range(ax.XLim);
    
    % The double headed arrow
    annotation('doublearrow', ax_pos(1) + [norm_start norm_end] * ax_pos(3), ax_pos(2) + [0.5 0.5] * ax_pos(4), 'Color', 'k');
    
    % The text annotation for each beam on duration
    text(mean([AP_time(onEdges(i)) AP_time(offEdges(i))]), 1.05, sprintf('%.2f s', beamOnDurations(i)), 'HorizontalAlignment', 'center');
end

% Adjust figure position if necessary
hh.Position = [163, 664, 1276, 438];
grid on, title(sprintf("Average of beam-on time: %1.2f s",avgDuration))


beamOnSignalSet = nan(length(locs)-1, max_length);
for iter1=1:length(locs)-1
    temp_data = beam_status(locs(iter1):locs(iter1+1));
    beamOnSignalSet(iter1, 1:length(temp_data)) = temp_data;
end
% saveas(hh, "????? ???? ?????????_2 (100ms ????).png")
%%
close all 
numPeriods = size(gatingSignalSet_normalized, 1);
hhh = figure,
hhh.Position = [204.7735849056604,227.4150943396226,1234.226415094339,874.5849056603772]
diffSetSt = []
diffSetEnd= []
for iter1=1:numPeriods
% for iter1=2:10
    tempData = gatingSignalSet_normalized(iter1, :);
    tempData_norm = (tempData - min(tempData(:))) / (max(tempData(:)) - min(tempData(:)));

    tempData_epid = epidSignalSet(iter1, :)

    originalIndicesPhase_temp = linspace(0, 100, length(tempData_norm));

    signalInterval = x_step_second; % ???? ???? ?????? ????????

    % Calculate the times at the 25th and 75th percentiles
    timeAt25 = signalInterval * (st_phase/100 * (length(tempData_norm) - 1));
    timeAt75 = signalInterval * (end_phase/100 * (length(tempData_norm) - 1));
    timeDurationPhase(iter1) = timeAt75-timeAt25;

    subplot(3,4,iter1), 
    plot(originalIndicesPhase_temp, tempData_norm, 'k--'), 
    grid on, 
    xlabel("Phase (%)"), 
    ylabel("Normalized amplitude"), 
    hold on;
    
    % Plot the beam status
    beamOnDataTemp = beamOnSignalSet(iter1, :);
    plot(originalIndicesPhase_temp, beamOnDataTemp, 'r--');

    tempData_epid_norm = (tempData_epid-min(tempData_epid))./(max(tempData_epid(:)-min(tempData_epid)));
    plot(originalIndicesPhase_temp, tempData_epid_norm, 'b--')
    
    % Add yellow box for beam on status
    beamOnIndices = find(beamOnDataTemp == 1);
    if ~isempty(beamOnIndices)
        yl = ylim; % Get the current y-axis limits
        % Ensure the patch has four x values and four corresponding y values
        patchX = [originalIndicesPhase_temp(beamOnIndices(1)), originalIndicesPhase_temp(beamOnIndices(end)), originalIndicesPhase_temp(beamOnIndices(end)), originalIndicesPhase_temp(beamOnIndices(1))];
        patchY = [yl(1), yl(1), yl(2), yl(2)];
        patch(patchX, patchY, 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Plot red star at the 25th and 75th percentiles
    % Ensure the indices used for stars are within the valid range
    starIndex25 = max(1, round(st_phase/100 * length(tempData_norm)));
    starIndex75 = min(length(tempData_norm), round(end_phase/100 * length(tempData_norm)));
    plot(originalIndicesPhase_temp(starIndex25), tempData_norm(starIndex25), 'r*');
    plot(originalIndicesPhase_temp(starIndex75), tempData_norm(starIndex75), 'r*');
    
    % Calculate and display the differences
    startBeamStatus = find(beamOnDataTemp == 1, 1, 'first');
    endBeamStatus = find(beamOnDataTemp == 1, 1, 'last');
    if ~isempty(startBeamStatus) && ~isempty(endBeamStatus)
        diffStart = (startBeamStatus-1).*x_step_second - timeAt25;
        diffEnd =  (endBeamStatus-1).*x_step_second - timeAt75;
        if(diffEnd < -0.1)
            % diffEnd = diffEnd +0.07;
        end
        % title(sprintf('Start Diff: %.2f s, End Diff: %.2f s', abs(diffStart), abs(diffEnd)));
        % diffSetSt = [diffSetSt; abs(diffStart)]
        % diffSetEnd= [diffSetEnd; abs(diffEnd)]

        title(sprintf('Start Diff: %.2f s, End Diff: %.2f s', diffStart, diffEnd));
        diffSetSt = [diffSetSt; diffStart]
        diffSetEnd= [diffSetEnd; diffEnd]
    else
        title('No Beam On Data');
    end

    
    hold off;
end

title_temp = sprintf('Averaged Start Diff.: %.2f s, Averaged End Diff.: %.2f s', mean(diffSetSt, 1), mean(diffSetEnd, 1))
sgtitle(title_temp)
% saveas(hhh, "????????? (100ms ????).png")
% sgtitle("hi")
% title("??????? ??: <0.1 s")

%%
% filenameFrameSum = "frame_sums.csv"
% FrameSumData = xlsread(filenameFrameSum )
% times = FrameSumData(:,1)
% frameSums = FrameSumData(:,2)
% frameSums_Binary = frameSums>2.35*10^8;
% figure, plot(times, frameSums)
% figure, plot(times, frameSums_Binary )
% 
% frameSums_Bianry_diff = diff(frameSums_Binary)
% 
% figure, plot(times(1:end-1), frameSums_Bianry_diff )
% 
% frameSums_Bianry_diff_1 = find(frameSums_Bianry_diff==1);
% frameSums_Bianry_diff_0 = find(frameSums_Bianry_diff==-1);
% frameSums_Bianry_diff_0(1) = []
% frameSums_Bianry_diff_1(end) = []
% 
% index_diff = frameSums_Bianry_diff_0-frameSums_Bianry_diff_1
% 
% frameRate = mean(diff(times))
% mean(index_diff.*frameRate)


