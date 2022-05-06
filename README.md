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





Usage
------------
Decide on a

License and Citation
------------




Definitions:
------------
Kisman-Masurel Original: 
| KM score | Bubbles per cardiac cycle | Percentage of cardiac cycles at rest with detectable bubbles at rest | Relative amplitude |
| --- | --- | --- | --- | 
| 0 | 0 | 0 | 0 |
| 1 | 1-2 | 1-10% | 0.100-0.250 |
| 2 | 3-8 | 10-50% | 0.250-0.450 |
| 3 | 9-P_y | 50-99% | 0.450-0.775 |
| 4 | P_y - (2* P_y) | 100% | 0.775-0.999 |

Kisman-Masurel Adapted: 
| KM score | Bubbles per cardiac cycle | Percentage of cardiac cycles at rest with detectable bubbles at rest | Relative amplitude |
| --- | --- | --- | --- | 
| 0 | 0 | 0 | 0 |
| 1 | 1-2 | 1-10% | 0.100-0.250 |
| 2 | 3-8 | 10-50% | 0.250-0.450 |
| 3 | 9-P_y | 50-99% | 0.450-0.775 |
| 4 | P_y - (2* P_y) | 100% | 0.775-0.999 |

P_y = {Cardiac Period (s)}/ {Average Bubble Length (s)}

Spencer Table:
| Grade | Description |
| --- | --- |
| Grade 0 | Complete lack of bubble signals|
| Grade 1 | Occasional bubble signal discernible with the cardiac motion signal, with majority of cardiac periods free of bubbles|
| Grade 2 | Many but less than half of the cardiac periods contain bubble signals, singu- larly or in groups|
| Grade 3 | All of the cardiac periods contain showers or single bubble signals, but not dominating or overriding the cardiac motion signals|
| Grade 4 | Maximum detectable bubble signal sounding continuously throughout systole and diastole of every cardiac period, and overriding the amplitude of the normal cardiac signals|

Kisman-Masurel to Spencer conversion Table: 
| Spencer | Kisman-Masurel Grades |
| --- | --- |
| 0 | 000 |
| 1 | 111 112 113 211 212 213 |
| 2 | 121 122 123 221 222 223 |
| 3 | 232 233 242 243 332 333 342 343 |
| 4 | 444 |
