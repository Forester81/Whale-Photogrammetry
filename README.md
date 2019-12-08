# Whale-Photogrammetry
Whale photogrammetry processing software for extracting scale morphometric attributes from baleen whale imagery

Updated: 12/08/19

Pre Processing Summary (Whale Measurement):
Navigate to this location to download the pre-processing software for scaling object and whale measurement:  https://oregonstate.app.box.com/s/tmiqgabvncoqvzi04a0imc2frs6307mj

Post Processing Summary (Measurement Analysis):
Whale_Calibration_Object_Measurement.m and Whale_Measurements.m are the main programs

use Whale_Calibration_Object_Measurement.m for compiling information about validation and training objects for calibration

use Whale_Measurements.m to measure morphometric properties of whales. 

The purpose of both programs is to standardize the way data is interpreted and output.  There is an analytical backend developed in R to consume the data products from both of these processes.

Calibration:

The program requires focal length (in mm) and pixel pitch (in mm) to estimate GSD. For DJI Phantom 3 Pro, Phantom 4 and Phantom 4 Pro calibration files have been developed. These calibration files are ONLY for video shot in 3840 x 2160 resolution only. All other formats and camera require a user provided Matlab calibration file OR the two parameters specified above. 

Custom Calibration File:

This program is designed to ingest a custom calibration file.  The calibration file must be created using Matlab's Camera Calibration application. To do this take at least 30 (recommend 50) images of the checkerboard pattern available in the app (or use Agisoft Lens' checkerboard). The checkerboard should be imaged from multiple angles and distances. Once the calibration coefficients have been created, output the results into the Matlab workspace with the name 'calibrationSession'.  It's important to note that the program is hardcoded to look for a workspace object named 'calibrationSession' that houses all of the calibration information. Any other name will throw an error. Save that workspace object into a .mat file with a name of your choosing.  Source this file by clicking the 'custom' option when running the program.   

Analysis:

'Whale_Quant_Analysis.R' is the number crunching backend to the program suite.  It ingests the outputs from the Whale_Measurements and Whale Calibration Object Measurement programs and produces a summary table. The summary table is morphometric attributes plus upper and lower confidence intervals, as well as coefficients of variation. The summary is broken out by each of the calibration methods described in the paper with Method 1 being the uncorrected 'as-measured' dataset. Summary results for the calibration object data are also output.  The output is a .csv places in the parent directory containing all of the image data.

Important:

'Whale_Quant_Funcs.R' contains all of the functions required by 'Whale_Quant_Analysis.R'.  Both files need to be located in the same folder so the 'Whale_Quant_Funcs.R' script can be sourced.

Tips:

Read prompts very carefully and click 'ok' when you are ready to proceed!  There is admittedly some minor glitchines s associated with the message boxes and the subsequent interactions, however, it's far superior to manual analaysis.  
