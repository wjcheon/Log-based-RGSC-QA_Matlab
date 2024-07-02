clc
clear
close all
%%
% filenameFrameSum = "\\10.10.37.115\Research\0.wjcheon\00_Research\07_GatingQA @ynkang\EPID video\2403\frame_sums2.csv";
filenameFrameSum = "\\10.10.37.115\Research\0.wjcheon\00_Research\07_GatingQA @ynkang\EPID video\trial_3\frame_sums_trail3.csv"
close all
FrameSumData = xlsread(filenameFrameSum );
times = FrameSumData(1:end,1);
frameSums = FrameSumData(1:end,2);
figure, plot(times, frameSums)

%Pre-processing
minVal = min(frameSums);
maxVal = max(frameSums);
% 최소값과 최대값을 뒤집기
frameSums = maxVal + minVal - frameSums;
% 정규화: 최대값이 1이 되도록 조정
frameSums = frameSums / max(frameSums);
figure, plot(times, frameSums)


% frameSums_Binary = frameSums>6.0 *10^6;
frameSums_Binary =frameSums> 0.8;
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





