# Synthetic-Post-Dive-Ultrasound-Audio-Generator

News
------------
2022/05/06: Repository created.
2022/05/~: Paper submitted??


Requirements
------------
MATLAB version 2020b or newer (needs audio data augmentation toolbox) 

Overview
------------
This repository contains the code described in Le et al. 2022, An open-source framework for synthetic post-dive Doppler ultrasound audio generation. 
Doppler ultrasound (DU) is recorded after divers surface to measure the amount of bubble (venous gas emboli or VGE) generation from decompression. These recordings are acquired of the subclavian vein or precordium (precordial). Precordial DU 

Large amounts of real-world data is difficult to acquire and do not provide ground truth data for what is VGE or cardiac signal. The algorithm presented here generates synthetic post-dive DU data using real-world clean baselines for precordial/subclavian vein recordings injected with isolated VGE audio signals recorded experimentally. A rule-based system is implemented for the injection of VGE into these data, following the Kisman-Masurel grading scale (Table 1,2). The KM scale can be converted into the Spencer scale (Table 3) which is also used by researchers using table 4. 

The code allows a user to generate any number of audio files using either Spencer or KM grading scales. Furthermore, the user can designate certain parameters such as length of audio file, sampling frequency. A method to place VGE only in quiet regions of a cardiac cycle is provided to create data with clearer separation between VGE and human signals. Finally, all code is modifiable such as the definition of KM scale. 

Example image of generated KM 222 data. 
![Example KM 222 image](https://github.com/dle4/Synthetic-Post-Dive-Ultrasound-Audio-Generator/blob/main/Main/Images/ExampleData.png)


Usage
------------

Two main functions are provided: 1) SyntheticCombinationDoppler_FullCardiac_Master_04_18.m and 2) SyntheticCombinationDoppler_PartialCardiac_Master_04_20.m. 
The fullcardiac script generates synthetic data where VGE are placed anywhere within detected cardiac cycles. PartialCardiac script only places VGE within the region between cardiac cycles using the FWHM of the peaks to define the placement of VGE. These scripts are directly used to generate data. To use, open the script in MATLAB and modify the section "User-defined parameters" with number of files to be generated (numfiles), desired audio file length in seconds (desired_length_sec), output data sampling frequency (Fs2), and code system (Spencer = 1, Kisman-Masurel = 2). Additionally, properly define the folder location of the baseline DU human and VGE audio files. Finally, specify folder and filename for generated data to be placed into and named. 

By modifying the source of the baseline human data, the user can choose between precordial or subclavian output data. 

Data is created following directory structure: 

![alt text](https://github.com/dle4/Synthetic-Post-Dive-Ultrasound-Audio-Generator/blob/main/Main/Images/ExampleDirectoryStructure.png | width = 100)

Example data is provided covering 6 cases ranging  subclavian/precordial, full/partial cardiac cycle, and spencer/KM scales. Partial cardiac cycle script is not used for subclavian data as the real-world data is often clean enough to easily differentiate background from VGE. 



License and Citation
------------




Definitions:
------------
Kisman-Masurel Original (Table 1): 
| KM score | Bubbles per cardiac cycle | Percentage of cardiac cycles at rest with detectable bubbles at rest | Relative amplitude |
| --- | --- | --- | --- | 
| 0 | 0 | 0 | No sound |
| 1 | 1-2 | 1-10% | Barely perceptible |
| 2 | 3-8 | 10-50% | Moderate amplitude|
| 3 | 9-P_y | 50-99% | Loud bubbles |
| 4 | P_y - (2* P_y) | 100% | Maximal|

Kisman-Masurel Adapted (Table 2): 
| KM score | Bubbles per cardiac cycle | Percentage of cardiac cycles at rest with detectable bubbles at rest | Relative amplitude |
| --- | --- | --- | --- | 
| 0 | 0 | 0 | 0 |
| 1 | 1-2 | 1-10% | 0.100-0.250 |
| 2 | 3-8 | 10-50% | 0.250-0.450 |
| 3 | 9-P_y | 50-99% | 0.450-0.775 |
| 4 | P_y - (2* P_y) | 100% | 0.775-0.999 |

P_y = {Cardiac Period (s)}/ {Average Bubble Length (s)}

Spencer Table (Table 3):
| Grade | Description |
| --- | --- |
| Grade 0 | Complete lack of bubble signals|
| Grade 1 | Occasional bubble signal discernible with the cardiac motion signal, with majority of cardiac periods free of bubbles|
| Grade 2 | Many but less than half of the cardiac periods contain bubble signals, singu- larly or in groups|
| Grade 3 | All of the cardiac periods contain showers or single bubble signals, but not dominating or overriding the cardiac motion signals|
| Grade 4 | Maximum detectable bubble signal sounding continuously throughout systole and diastole of every cardiac period, and overriding the amplitude of the normal cardiac signals|

Kisman-Masurel to Spencer conversion Table (table 4): 
| Spencer | Kisman-Masurel Grades |
| --- | --- |
| 0 | 000 |
| 1 | 111 112 113 211 212 213 |
| 2 | 121 122 123 221 222 223 |
| 3 | 232 233 242 243 332 333 342 343 |
| 4 | 444 |
