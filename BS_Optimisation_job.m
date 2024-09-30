%-----------------------------------------------------------------------
% Job saved on 09-Jul-2024 17:10:52 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (12.6)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

% Updated 24/09/2024
% by Shokoufeh Golshani


toolboxFileLoc = [fullfile(spm('Dir'), 'toolbox','FmpOptBS') filesep];

% =========================================================================
% Add the folder containing the fieldmap or fieldmap gradients (*.nii)
% =========================================================================
datapath = uigetdir([],'Select the fieldmaps Directory'); 
files = dir(fullfile(datapath, '*.nii'));

Fieldmap_files = {files(~[files.isdir]).name};
charArray = char(Fieldmap_files);

inpt = [fullfile(toolboxFileLoc,'Fieldmap_grads') filesep];

n_fmap = size(charArray,1);
fmap = [];
for i = 1:n_fmap
    fmap{i,1} = [inpt charArray(i,:)];
end

matlabbatch{1}.spm.tools.FmpOptBS.inputfiles.fieldmaps = fmap;


% =========================================================================
% Add the folder containing the BrainMask template (*.nii)
% =========================================================================
datapath = uigetdir([],'Select the template Directory'); 
files = dir(fullfile(datapath, '*.nii'));

Temp_files = {files(~[files.isdir]).name};
charArray = char(Temp_files);

inpt = [fullfile(toolboxFileLoc,'Template') filesep];

n_tmpl = size(charArray,1);
tmpl = [];
for i = 1:n_tmpl
    tmpl{i,1} = [inpt charArray(i,:)];
end

matlabbatch{1}.spm.tools.FmpOptBS.inputfiles.template = tmpl;


% =========================================================================
% Add the folder containing the desired brain ROIs (*.nii)
% =========================================================================
datapath = uigetdir([],'Select the ROIs Directory'); 
files = dir(fullfile(datapath, '*.nii'));

ROI_files = {files(~[files.isdir]).name};
charArray = char(ROI_files);

inpt = [fullfile(toolboxFileLoc,'ROIs') filesep];

n_ROI = size(charArray,1);
Myroi = [];
for i = 1:n_ROI
    Myroi{i,1} = [inpt charArray(i,:)];
end

matlabbatch{1}.spm.tools.FmpOptBS.inputfiles.rois = Myroi;

