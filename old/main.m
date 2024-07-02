clc
clear
close all
fclose all
%%

filename = '18815802_1.vxp';
fileID = fopen(filename, 'r');

% 빈 구조체 초기화
dataStruct = struct('amplitude', [], 'phase', [], 'timestamp', [], 'validflag', [], 'ttlin', [], 'mark', [], 'ttlout', []);

% 파일에서 한 줄씩 읽기
line = fgetl(fileID);
line = fgetl(fileID);
line = fgetl(fileID);
line = fgetl(fileID);
line = fgetl(fileID);
line = fgetl(fileID);
line = fgetl(fileID);
line = fgetl(fileID);
numStr = regexp(line, '\d+', 'match');
samplePerSecond = str2double(numStr);

line = fgetl(fileID);
numStr = regexp(line, '\d+\.?\d*', 'match');
scaleFactor = str2double(numStr);

line = fgetl(fileID);
while ischar(line)
    line = fgetl(fileID);
    if line == -1
        % 파일의 끝
        break;
    end
    % 쉼표로 구분하여 데이터 파싱
    a=0;

    [C,matches] = strsplit(line, ",");

    % 구조체에 데이터 추가
    if ~isempty(C)
        dataStruct.amplitude(end+1) = str2double(C{1});
        dataStruct.phase(end+1) = str2double(C{2});
        dataStruct.timestamp(end+1) = str2double(C{3});
        dataStruct.validflag(end+1) = str2double(C{4});
        dataStruct.ttlin(end+1) = str2double(C{5});
        dataStruct.mark(end+1) = str2double(C{6});
        % dataStruct.ttlout(end+1) = str2double(C{7});
    end

end

fclose(fileID);
disp(dataStruct);
%%
close all
% figure, 
% subplot(2,3,1), plot(dataStruct.amplitude)
% subplot(2,3,2), plot(dataStruct.phase)
% subplot(2,3,3), plot(dataStruct.timestamp)
% subplot(2,3,4), plot(dataStruct.validflag)
% subplot(2,3,5), plot(dataStruct.ttlin)
% subplot(2,3,6), plot(dataStruct.mark)

data_norm = (dataStruct.amplitude)./max(dataStruct.amplitude(:))
data_norm_smooth = smoothdata(data_norm, 5)

[pks,locs,w,p] = findpeaks(data_norm, 'MinPeakDistance',60, 'MinPeakProminence', 0.01)
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
for iter1=1:length(locs)-1
    temp_data = data_norm_smooth(locs(iter1):locs(iter1+1));
    gatingSignalSet(iter1, 1:length(temp_data)) = temp_data;
end

gatingSignalSetAvg = nanmean(gatingSignalSet, 1)

figure, plot(gatingSignalSetAvg)
%%
% 신호를 10개의 phase로 나눔
gatingSignalSetAvg = (gatingSignalSetAvg - min(gatingSignalSetAvg)) / (max(gatingSignalSetAvg) - min(gatingSignalSetAvg));
gatingSignalSetAvg = gatingSignalSetAvg.*100;
numPhases = 10;
phaseLength = length(gatingSignalSetAvg) / numPhases;
% 신호를 10개의 phase로 나눔
XData = 1:length(gatingSignalSetAvg)

% 각 phase 처리
phaseSet = []
ampSet = []
for phase = 1:numPhases
    % 현재 phase의 신호 부분과 해당 시간 벡터 추출
    startIndex = round((phase-1)*phaseLength) + 1;
    endIndex = round(phase*phaseLength);
    currentPhaseSignal = gatingSignalSetAvg(startIndex:endIndex);
    currentPhaseX = XData(startIndex:endIndex); % 시간 벡터도 추출
    
    phaseSet(phase)= startIndex;
    ampSet(phase) = gatingSignalSetAvg(startIndex);

    meanValue = mean(currentPhaseSignal);
    fprintf('Phase %d: Mean Value = %f\n', phase, meanValue);
    
    % 필요한 경우 현재 phase의 시간 정보와 신호 정보를 이용한 추가 처리
    % 예: plot 그리기
    % plot(currentPhaseTime, currentPhaseSignal);
    % hold on; % 모든 phase를 같은 그래프에 표시
end
% hold off; % 모든 phase의 그래프가 그려진 후 hold를 해제
%% 
close all
figure, plot(gatingSignalSetAvg), hold on, grid on 
plot(phaseSet, ampSet, 'r*')
xlabel('Pixel'), ylabel('Normalized amplitude')
