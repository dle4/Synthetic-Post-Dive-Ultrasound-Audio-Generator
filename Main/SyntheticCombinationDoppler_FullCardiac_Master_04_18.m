%% Version 5 of the synthetic data code
% written by David Le 03/25/2022
% edited 04/18/2022
% This script functions by allowing a user to generate synthetic Doppler
% data containing Venous Gas Emboli in Precordial Doppler Ultrasound. 
% The VGE are placed into pre-existing baseline cardiac Doppler audio based 
% on the Kisman-Masurel scale and subsequently sorted into Spencer grades 
% according to a conversion chart. 

% no additional scripts are needed. Can change the code to be parallized. 
% using parfor instead of for on line 154

% To use the script, specify user-defined variables such as number of files
% per class, desired audio length (in seconds), and final sampling
% frequency. Also, Kisman-Masurel or Spencer style data output. 

% File will be generated with name "UserDefinedName_Class_index.wav" and
% placed into relevant class folders. 

% Three separate folders groupings will be generated, 

%% housekeeping 
clear; close all; clc; 

%% User-defined parameters
numfiles = 1; %number of synthetic Doppler files per class to be made
desired_length_sec = 10; %seconds per audio file
Fs2 = 8000; %resampled frequency for output data (44100 and 8000 are good choices)
codesystem = 2; % 1 for Spencer and 2 for Kisman-Masurel

% Specify folders where baseline data is located
% baseline_human_dir = "D:\Projects\Doppler Project\Data\Simulated data\Rawdata\O'Dive dataset - only pre-dive"; 
baseline_human_dir = 'D:\Projects\Doppler Project\Data\Simulated data\Rawdata\2021dopplercardiac072921';
bubble_dir = 'D:\Projects\Doppler Project\Data\Simulated data\Rawdata\SimulatedBubbles_Sequoia';

% where to save augmented data
savefolder_all = 'E:\Projects\Doppler Project\Data\Simulated data\Synthetic Doppler Data\TestBubbles_Spencer_10s_2022_04_18\';
savefolder_cardiac = [savefolder_all 'DopplerSynthCardiac\'];
savefolder_bubbles = [savefolder_all 'DopplerSynthBubbles\'];
savefolder_combined = [savefolder_all 'DopplerSynthCombined\'];
savefilebasename = 'syntheticDopplerAudioCombined_'; %what the file will be named

try
    mkdir(savefolder_cardiac);
    mkdir(savefolder_bubbles);
    mkdir(savefolder_combined);
end

%% Definitions for Kisman Masurel code (grouped for Spencer)
switch codesystem
    case 1 % Spencer code
        sf = {};
        sf{1} = [0 0 0];
        sf{2} = [1 1 1; 1 1 2; 1 1 3; 2 1 1; 2 1 2; 2 1 3];
        sf{3} = [1 2 1; 1 2 2; 1 2 3; 2 2 1; 2 2 2; 2 2 3];
        sf{4} = [2 3 2; 2 3 3; 2 4 2; 2 4 3; 3 3 2; 3 3 3; 3 4 2; 3 4 3];
        sf{5} = [4 4 4];
        all_classes = {'1','2','3','4','5'}; 
    case 2 % Kisman Masurel code, but only the ones that are possible according to Bennett and Elliot
        sf = {};
        sf{1} = [0 0 0];
        sf{2} = [1 1 1; 1 1 2; 1 1 3; 2 1 1; 2 1 2; 2 1 3];
        sf{3} = [1 2 1; 1 2 2; 1 2 3; 2 2 1; 2 2 2; 2 2 3];
        sf{4} = [2 3 2; 2 3 3; 2 4 2; 2 4 3; 3 3 2; 3 3 3; 3 4 2; 3 4 3];
        sf{5} = [4 4 4];

        all_sf = cat(1,sf{:});
        sf =num2cell(all_sf,2)';
        for i = 1:length(all_sf)
            name = num2str(all_sf(i,:));
            name(name == ' ') = [];
            all_classes{i} = name;
        end
    otherwise
        disp('Choose 1 or 2')
end
%% load bubble only data and extract single bubbles from the files
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
                individual_bbl2 = resample(individual_bbl, Fs2,Fs); %resample bubbles to desired Fs2
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


%% Load human Doppler data into matlab and preprocess
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




%% Main loop of code that generates data
% change from "for" to "parfor" if parallelization is desired
parfor f = 1:length(sf)
    class_name = all_classes{f};
    fpname1 = [savefolder_cardiac, class_name];
    fpname2 = [savefolder_bubbles, class_name];
    fpname3 = [savefolder_combined, class_name];

    mkdir(fpname1);
    mkdir(fpname2);
    mkdir(fpname3);


%     direc_check = dir(fpname3);
    while (length(dir(fpname3))-2) < numfiles
        try
            %% Randomly sample an audio signal from cardiac data and augment without noise
            randaudio = audioAll{randi(length(audioAll))};
            desired_length_samp = desired_length_sec*Fs2;
            buffered_length = desired_length_samp*1.5;
            start_loc = randi(length(randaudio)-buffered_length);
            randsample = randaudio(start_loc:start_loc+buffered_length-1);
            randsample = randsample./max(randsample);
            %%
            % augment the audio
            augmenter = audioDataAugmenter( ...
                "AugmentationMode","sequential", ...
                "NumAugmentations",1, ...
                ...
                "TimeStretchProbability",0.5, ...
                "SpeedupFactorRange", [0.8,1.3], ...
                ...
                "PitchShiftProbability",0.5, ...
                ...
                "VolumeControlProbability",0, ...
                "VolumeGainRange",[-10,10], ...
                ...
                "AddNoiseProbability",0, ...
                "SNRRange", [0,2], ...
                ...
                "TimeShiftProbability",0.5, ...
                "TimeShiftRange", [-desired_length_sec, desired_length_sec]);

            data = augment(augmenter,randsample,Fs2);
            rand_sample_aug = data.Audio{1};

            %% Determine sections of the cardiac signal that are available to place bubbles
            curr_seg = rand_sample_aug;
            close
            y = curr_seg./max(curr_seg); %normalize signal

            y2 = abs(hilbert(y)); %envelope detect
            y2 = movmean(-1*y2,500); %average filter and flip the signal
            y2 = y2-min(y2); %shift the signal to positive
            y2 = y2./max(y2); %normalize

            %apply hrd to determine minpeakdistance and npeaks
            Nmin = (60*8000)/120; % minimum sample index (searching window for both)
            Nmax = (60*8000)/50;  % maximum sample index
            [autocor,lags] = xcorr(y2, 'unbiased');
            autocor(1:length(y2)-1) = [];  % removing the unnecessary part
            lags(1:length(y2)-1) = [];
            lagfinalInstHR = lags;
            lagfinalInstHR(Nmax:end) = [];
            lagfinalInstHR(1:Nmin) = [];
            autocor2 = autocor;
            autocor2(Nmax:end) = [];
            autocor2(1:Nmin) = [];
            [pks1, locs1] = findpeaks(autocor2);
            [maxVal, maxLoc] = max(pks1);
            minpeakdist = lagfinalInstHR(locs1(maxLoc));
            npeaks = round(length(y2)/minpeakdist);
            heartrate = (60*Fs2) / lagfinalInstHR(locs1(maxLoc));


            %perform peak detection to determine where heartbeats occur
            [pks2, locs2] = findpeaks(y2,"NPeaks",npeaks, "MinPeakDistance",minpeakdist*.9); % detect troughs in the inverted cardiac signal

            window_indexes = [locs2(1:end-1) locs2(2:end)];%since we have detected the troughs of the signal, we will shift all windows by 50% to capture the peaks
            window_shifts = (locs2(2:end)-locs2(1:end-1))/2;
            window_indexes2 =window_indexes+ window_shifts;
            indices = find(sum(window_indexes2>length(curr_seg),2),1); % remove any indices that are outside the actual length of the signal
            window_indexes2(indices,:) = [];


            avg_bbl_lgnth_sec = avg_bbl_lngth/Fs2;
            audio_length = length(rand_sample_aug);
            % conversion between KM and Spencer

            km_choice = sf{f};
            km = km_choice(randi([1 size(km_choice,1)]),:);

            selected_cycles = [];
            if km ~= [0,0,0] % create the bubble_audio here
                bubble_audio = zeros(audio_length,1)+(rand(audio_length,1)-0.5)*1e-4;
                % determine percentage of cardiac cycles at rest with detectable bubbles
                if km(2)==1
                    per_cardiac_cycles = round(size(window_indexes2,1)*randi([1,10])/100);
                    if  per_cardiac_cycles == 0
                        per_cardiac_cycles = 1;
                    end
                elseif km(2)==2
                    per_cardiac_cycles = round(size(window_indexes2,1)*randi([10,50])/100);
                elseif km(2)==3
                    per_cardiac_cycles = round(size(window_indexes2,1)*randi([50,99])/100);
                elseif km(2)==4
                    per_cardiac_cycles = size(window_indexes2,1);
                end


                selected_cycles = datasample(window_indexes2,per_cardiac_cycles,'Replace',false); %randomly sample a certain number of cardiac cycles from the cardiac signal (in win_idx)

                % amplitude
                if km(3)==1
                    amplitude = randi([100,250])/1000;
                elseif km(3)==2
                    amplitude = randi([250,450])/1000;
                elseif km(3)==3
                    amplitude = randi([450,775])/1000;
                elseif km(3)==4
                    amplitude = randi([775,999])/1000;
                end
                for j = 1:size(selected_cycles,1) % iterate through all the windows that were selected and define the number of bbls that will be placed in them.
                    cycle_curr_length = selected_cycles(j,2)-selected_cycles(j,1); % end of window - start of window to define length
                    cardiac_period = cycle_curr_length/Fs2; %convert number of samples to time in (s)

                    if km(1)==1 %using km1, define number of bubbles per cardiac cycle
                        num_bbls = randi([1,2]);
                    elseif km(1)==2
                        num_bbls = randi([3,8]);
                    elseif km(1)==3
                        if round(cardiac_period/avg_bbl_lgnth_sec) < 9
                            num_bbls = 9;
                        else
                            num_bbls = randi([9,round(cardiac_period/avg_bbl_lgnth_sec)]);
                        end
                    elseif km(1)==4
                        if round(cardiac_period/avg_bbl_lgnth_sec) < 9 % this ensures that greater than 9 bubbles are placed into the section when the cardiac period is too short for 9
                            num_bbls = randi([12,18]);
                        else
                            num_bbls = randi([round(cardiac_period/avg_bbl_lgnth_sec),round(cardiac_period/avg_bbl_lgnth_sec)*2]); %this calculates # of bubbles needed for continuous sound (and up to double with overlap)

                        end
                    end


                    sampled_bbls = datasample(bubble_seg_resample,num_bbls); %randomly sample num_bbls from the all individual bubbles found in "signal_seg_resample"

                    %augment the individual bubbles
                    num_augments = 1;
                    augmenter_bubble = audioDataAugmenter( ...
                        "AugmentationMode","sequential", ...
                        "NumAugmentations",num_augments, ...
                        "TimeStretchProbability",0.5, ...
                        "SpeedupFactorRange", [0.8, 2], ...
                        "PitchShiftProbability",.5, ...
                        "VolumeControlProbability",0, ...
                        "AddNoiseProbability",0, ...
                        "TimeShiftProbability",0);
                    bubble_audio_singlecycle = zeros(cycle_curr_length,1)+(rand(cycle_curr_length,1)-0.5)*1e-4; %allocate the single cardiac cycle bubble audio array

                    for i = 1:length(sampled_bbls) %iterate through each bubble, augmenting them and then placing into the signal.
                        curr_bbl = sampled_bbls{i};
                        if length(curr_bbl)< 1323  % this catches errors for the augmentation. Signals less than 1323 cannot be fed into it.
                            fixed_length = zeros([1323-length(curr_bbl),1]);
                            nbbl = [curr_bbl;fixed_length];
                            curr_bbl = nbbl;
                        end
                        norm_bbl = curr_bbl./max(curr_bbl);
                        % calculate amplitude of individual bubbles based on km3 (slightly
                        % random each bubble.

                        rand_amp = amplitude*(1+randi([-150 150])/1000); %introduce +- 15% variation to amplitude
                        amped_sig = rand_amp*norm_bbl;
                        data = augment(augmenter_bubble,amped_sig,Fs2);
                        data = augment(augmenter_bubble,amped_sig,Fs2);
                        aug_sig = data.Audio{1};

                        if length(aug_sig) < cycle_curr_length
                            bubble_drop = randi([1 cycle_curr_length-length(aug_sig)]);  %calculate where inthe cardiac cycle to place the bubble (randomly)
                        else
                            bubble_drop = 1;
                        end
                        pre_sig = zeros([bubble_drop,1]); % make sure the arrays are same size before adding new bubble to singlecycle bubble audio
                        post_sig = zeros([cycle_curr_length-bubble_drop-length(aug_sig),1]);
                        new_sig = [pre_sig; aug_sig; post_sig];
                        bubble_audio_singlecycle = bubble_audio_singlecycle+new_sig;

                    end

                    bubble_audio(selected_cycles(j,1):selected_cycles(j,1)+length(bubble_audio_singlecycle)-1) = bubble_audio_singlecycle;  %replace empty regions of bubble_audio with bubbles
                end

            else
                bubble_audio = zeros(audio_length,1)+(rand(audio_length,1)-0.5)*1e-4;
            end

            cardiac_audio = y; %make sure to use the normalized cardiac data
            combined_audio = bubble_audio+cardiac_audio;
            combined_audio = combined_audio./max(abs(combined_audio));  %normalize the output data to avoid clipping
            cropping = (audio_length-desired_length_samp)/2;

            bubble_audio = bubble_audio(cropping:desired_length_samp+cropping-1); % crop all the data to desired length by saving the center of the signal
            combined_audio = combined_audio(cropping:desired_length_samp+cropping-1);
            cardiac_audio = cardiac_audio(cropping:desired_length_samp+cropping-1);


            count = length(dir(fpname3))-1;
            savefilename_cardiac= [fpname1, '\',savefilebasename,class_name,'_'  num2str(count), '.wav'];
            savefilename_bubbles= [fpname2, '\',savefilebasename,class_name,'_'  num2str(count), '.wav'];
            savefilename_combined= [fpname3, '\',savefilebasename,class_name,'_'  num2str(count), '.wav'];
            audiowrite(savefilename_cardiac,cardiac_audio,Fs2)
            audiowrite(savefilename_bubbles,bubble_audio,Fs2)
            audiowrite(savefilename_combined,combined_audio,Fs2)

        end
    end


end

%% helper functions
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