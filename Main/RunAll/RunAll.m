%% housekeeping 
clear; close all; clc; 

%% User-defined parameters
numfiles = 1000; %number of synthetic Doppler files per class to be made
desired_length_sec = 10; %seconds per audio file
Fs2 = 8000; %resampled frequency for output data (44100 and 8000 are good choices)

%% Specify folders where baseline data is located
baseline_human_dir_precordial = 'D:\Projects\Doppler Project\Data\Simulated data\Rawdata\2021dopplercardiac072921';
baseline_human_dir_subclavian = "D:\Projects\Doppler Project\Data\Simulated data\Rawdata\O'Dive dataset - only pre-dive"; 
bubble_dir = 'D:\Projects\Doppler Project\Data\Simulated data\Rawdata\SimulatedBubbles_Sequoia';

[bubble_seg_resample,avg_bbl_lngth] = load_bbl_data(bubble_dir,Fs2); %load bubble data once
%% where to save augmented data
savefolder_all = 'E:\Projects\Doppler Project\Data\Simulated data\Synthetic Doppler Data\SyntheticDU_examples\';


%% Precordial
audioAll = load_humanDU_audio(baseline_human_dir_precordial,Fs2);
%
% case 1
savefolder_case1 = ['Spencer_Precordial_FullCardiacCycle\'];
savefilebasename = 'Spencer_Precordial_FullCardiacCycle_';
[savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case1);
codesystem = 1; % 1 for Spencer and 2 for Kisman-Masurel
SyntheticCombinationDoppler_FullCardiac_RunMulti_04_18; 
disp("case 1 done")
%
% case 2
savefolder_case2 = ['Spencer_Precordial_PartialCardiacCycle\'];
savefilebasename = 'Spencer_Precordial_PartialCardiacCycle_';
[savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case2);
codesystem = 1; % 1 for Spencer and 2 for Kisman-Masurel
SyntheticCombinationDoppler_PartialCardiac_RunMulti_05_06; 
disp("case 2 done")
% case 3
savefolder_case3 = ['KismanMasurel_Precordial_FullCardiacCycle\'];
savefilebasename = 'KM_Precordial_FullCardiacCycle_';
[savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case3);
codesystem = 2; % 1 for Spencer and 2 for Kisman-Masurel
SyntheticCombinationDoppler_FullCardiac_RunMulti_04_18; 
disp("case 3 done")
% case 4
savefolder_case4 = ['KismanMasurel_Precordial_PartialCardiacCycle\'];
savefilebasename = 'KM_Precordial_PartialCardiacCycle_';
[savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case4);
codesystem = 2; % 1 for Spencer and 2 for Kisman-Masurel
SyntheticCombinationDoppler_PartialCardiac_RunMulti_05_06; 
disp("case 4 done")

%% Subclavian
audioAll = load_humanDU_audio(baseline_human_dir_subclavian,Fs2); 
%
% case 5
savefolder_case5 = ['Spencer_Subclavian_FullCardiacCycle\'];
savefilebasename = 'Spencer_Subclavian_FullCardiacCycle_';
[savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case5);
codesystem = 1; % 1 for Spencer and 2 for Kisman-Masurel
SyntheticCombinationDoppler_FullCardiac_RunMulti_04_18; 
disp("case 5 done")
% case 6
savefolder_case6 = ['KismanMasurel_Subclavian_FullCardiacCycle\'];
savefilebasename = 'KM_Subclavian_FullCardiacCycle_';
[savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case6);
codesystem = 2; % 1 for Spencer and 2 for Kisman-Masurel
SyntheticCombinationDoppler_FullCardiac_RunMulti_04_18; 
disp("case 6 done")




%% Functions

function [savefolder_cardiac, savefolder_bubbles, savefolder_combined] = makesavefolder(savefolder_all,savefolder_case)

savefolder_cardiac = [savefolder_all savefolder_case 'DopplerSynthCardiac\'];
savefolder_bubbles = [savefolder_all savefolder_case 'DopplerSynthBubbles\'];
savefolder_combined = [savefolder_all savefolder_case 'DopplerSynthCombined\'];


try
    mkdir(savefolder_cardiac);
    mkdir(savefolder_bubbles);
    mkdir(savefolder_combined);
end

end


function [bubble_seg_resample,avg_bbl_lngth] = load_bbl_data(bubble_dir,Fs2)

% load bubble only data and extract single bubbles from the files
% load bubbles from a .mat file generated using SegmentBubbleSim_allFiles.m
direc_bbl = dir(bubble_dir);
counter = 1;
buffer = 0.01; % amount of buffer space around each detected bubble to be saved in (seconds)
bubble_seg = {};
bubble_seg_resample = {};
for nmx = 3:size(direc_bbl,1)  %this is the "set"
    filepath1 = [direc_bbl(nmx).folder '/' direc_bbl(nmx).name];
    [audio, Fs] = audioread(filepath1);
    audio1 = audio-mean(audio); %recenter data
    audio2 = movmax(audio1,round(Fs/20));
    audio3 = movmedian(audio2,round(Fs/80));
    audio4 = audio3/max(audio3);

    [pks, locs] = findpeaks(audio4,Fs,'MinPeakDistance',.1,'MinPeakHeight',.1, 'MinPeakProminence',.01);
    window_idx = [];
    bufferFs = buffer*Fs;
    for i = 1:size(pks,1)
        locPeak = round(locs(i)*Fs);
        vPeak = pks(i);
        idx_start = find(audio4(1:locPeak) < 0.2*vPeak, 1, 'last');
        idx_end = locPeak + find(audio4(locPeak:end) < 0.2*vPeak, 1, 'first');

        if idx_end-idx_start < 0.5*Fs  %no bubbling allowed to be longer than 0.5 seconds (maybe too long?)
            window_idx(i,1) = idx_start;
            window_idx(i,2) = idx_end;
            if idx_start-bufferFs > 0
                individual_bbl = audio(idx_start-bufferFs: idx_end+bufferFs);
                individual_bbl2 = resample(individual_bbl, Fs2,Fs);
                bubble_seg{counter} = individual_bbl;
                bubble_seg_resample{counter} = individual_bbl2;
                ind1(counter)= length(individual_bbl2);
                counter = counter+1;
            end
        end
    end
end

avg_bbl_lngth = mean(ind1);

disp('Finished Bubble Segmentation');
end

% Load human Doppler data into matlab and preprocess
function audioAll = load_humanDU_audio(baseline_human_dir,Fs2)

direc_baseline = dir(baseline_human_dir);
allnames = {};
allnames = checkdir(allnames,direc_baseline);
audioAll = {}; %preallocate the variable
for nmx = 1:size(allnames,2)
    filepath = allnames{nmx};
    [audioIn, Fs] = audioread(filepath);

    v{1} = audioIn(:,1);    % remove channels that are low energy
    idx = 1;
    if size(audioIn,2) == 2
        v{2} = audioIn(:,2);
        u1 = mean(abs(hilbert(v{1})));
        u2 = mean(abs(hilbert(v{2})));

        if u1 > u2
            idx = 1;
        elseif u2 > u1
            idx = 2;
        end
    end
    audioraw{nmx} = audioIn;
    audio1 = v{idx};

    curr_audio = resample(audio1, Fs2, Fs); %resample data to desired sampling frequency
    audioAll{nmx} = curr_audio; %save a single channel of data
    lngths(nmx) = length(audioIn);
end

disp('Finished Human Data Preprocess');
end

function [allnames] = checkdir(allnames,direc)
    %this function recursively searches through all subfolders and collects
    %files
    for i = 3:length(direc)
        count = length(allnames)+1; 
        if direc(i).isdir == 0
            allnames{count} = [direc(i).folder ,'/', direc(i).name];
        elseif direc(i).isdir==1
            direc_new = dir([direc(i).folder '/' direc(i).name]);
            [allnames]=checkdir(allnames,direc_new);
        end
            
    end

end