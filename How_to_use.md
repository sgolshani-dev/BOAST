How to Use:

-- Ensure SPM12 is installed (download from https://www.fil.ion.ucl.ac.uk/spm/software/download/).
   
-- Copy the FmpOptBS toolbox into the ~/SPM/toolbox directory.

-- Add the toolbox directory and all subfolders to the MATLAB path.

-- Copy your fieldmap gradients into the folder named Fieldmap_grads.
If you haven't calculated the gradients, copy the field map itself, and the code will compute the gradients automatically.

-- Copy your brain mask into the folder named Template.

-- Copy your ROIs into the folder named ROIs.

-- Launch SPM for fMRI analysis.

-- Open the batch interface and select BS_Optimisation_job.
Choose the appropriate input files as prompted.

-- Press the play button to start the process.

If you need to analyze the BS.matrix and your selected ROIs, add a pause at line 209 in the epi_opt_param_TB script.


Modifying Parameters:

The following scripts contain fixed protocol parameters that you can adjust as needed:

- SetDefaultEPIparam: Fixed parameters for the EPI protocol.
- SetDefaultScannerParam: Fixed parameters for the scanner settings.
- SetDefaultSimulationParam: Fixed parameters for the simulation.


