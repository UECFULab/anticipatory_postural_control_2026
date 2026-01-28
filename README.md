### Overview 
Data  & Program files for<br/> 
**Anticipatory Postural Control Emerges from a Predictive and Optimized Strategy for Movement Preparation**<br/>
Authors: Tetsuro Funato, Miho Ogawa, Akira Konosu, Dai Yanagihara

### System Requirements 
The analysis and simulation software was developed and tested using MATLAB R2024b. To run the program code successfully, the following MATLAB toolboxes are required: MATLAB Model Predictive Control Toolbox and Statistics and Machine Learning Toolbox.

## Data and program files for Experimental studies (in Analysis Folder)
### Structure of the measured data files

<ins>**Analysis/data/MotionData_subject[subject number].xlsx**</ins>

* **comX,** **angAnkle,** **angHip,** **angForcePlate,** **posAnkleX:** Time series of motion data. Rows correspond to time, and columns correspond to trials (60 trials in total). Each column shows the time series of motion measured by a motion capture system for each trial. The time series range in the data is 5 seconds before and after the start of the floor tilting. The measurement frequency is 300 Hz. "comX": horizontal COM displacement (mm), "angAnkle": ankle angle (rad), "angHip": hip angle (rad), "angForcePlate": force plate angle (rad), and "posAnkleX": horizontal ankle displacement (mm). 
* **tiltEndIndx,** **Soundstart,** **Soundend:** Time of end of floor tilting, start of cue, and stop of cue. Each column corresponds to trial (60 trials in total). Each value is represented as the index number of the measured data. "tiltEndIndx" is measured with 300 Hz, "Soundstart" and "Soundend" are measured with 1000 Hz.The actual time (seconds) is calculated by dividing the data index by 300 or 1000.

<ins>**Analysis/data/EMGData_subject[subject number].xlsx**</ins>

* **RTA,** **LTA,** **RGAS,** **LGAS:** Time series of Electromyography (EMG) data. Rows correspond to time, and columns correspond to trials (60 trials in total). Each column shows the time series of EMG for each trial. The time series range in the data is 5 seconds before and after the start of the floor tilting. The measurement frequency is 300 Hz.  Each value represents %MVC, which is the measured EMG divided by maximum muscle activity. "RTA": right tibialis anterior (TA), "LTA": left TA, "RGAS": right gastrocnemius (GC), and "LGAS": left GC.

### Structure of the program files (Matlab files)
<ins>**Analysis/analysisMotion.m**</ins>

The program loads measured motion data and analyzes the time series of each subject's movements and differences in movements between cue and no cue tasks.

* Input data: **MotionData_subject[subject number].xlsx**
* Output data:
  * **tResult{COM/Hip/Ank}.csv:** Results of t-tests to determine significant differences in COM/hip/ankle changes in the CS-FS interval between cue and no cue tasks.
* Output figures:
  * **{COM/Hip/Ankle}_subject{subject_number}.emf:** Time series of COM/Hip/Ankle.
  * **COMAnkleHipdiff.emf:** The amount of change in COM/hip/ankle values for all subjects between CS and FS.

<ins>**Analysis/analysisEMG.m**</ins>

The program loads measured EMG data and analyzes the time series of each subject's EMGs and differences in EMGs between cue and no cue tasks.

* Input data: **EMGData_subject[subject number].xlsx**
* Output data:
  * **tResult{TA/GC}.csv:** Results of t-tests to determine significant differences in TA/GC activities in the CS-FS interval between cue and no cue tasks.
* Output figures:
  * **{TA/GC}_subject{subject_number}.emf:** Time series of TA/GCe.
  * **cTAGC_ALL.emf:** Coefficient for linear regression of TA/GC activities in the CS-FS interval for all subjects.


## Program files for Simulation studies (in Simulation Folder)
### Structure of the simulation program
<ins>**stand_ms.m**</ins>

Main program.

* Output files:
  * **simResults/simResultData_{simulation number}.csv:** Simulation results (Time series of COM, Joint angles, Muscle activation, Floor angles, etc.)
  * **simResults/simInfo_{simulation number}.csv:** Properties used for simulation (Seed of the random number, COM Reference, Weight for the evaluation function) 

<ins>**diffEqLink2_ms.m**</ins>

Differential equation of the system

<ins>**odeEulerModel.m**</ins>

Internal model used for MPC

<ins>**standOutFunc_ms.m**</ins>

Output function used for MPC 

### Structure of the result display program
<ins>**showResultFromData.m**</ins>

Load simulation results in "simResults" folder and display results.
