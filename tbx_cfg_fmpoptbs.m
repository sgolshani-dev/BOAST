function fmpoptbs = tbx_cfg_fmpoptbs

% =========================================================================
% MATLABBATCH Configuration file for toolbox 'tbx_cfg_FmpOptBS'
%_______________________________________________________________________
% Copyright (C) 2015-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 

% $Id: tbx_cfg_FmpOptBS.m 89 2015-12-16 11:03:00Z steffen $
% =========================================================================

% Updated 23/09/2024
% by Shokoufeh Golshani


% Adding the toolbox folder
if ~isdeployed
    addpath(fullfile(spm('Dir'),'toolbox','FmpOptBS')); 
end

%==========================================================================
% Default values that are common to all tbx_cfg_FmpOptBS jobs
%==========================================================================

% -------------------------------------------------------------------------
% field maps
% -------------------------------------------------------------------------
fieldmaps         = cfg_files;
fieldmaps.tag     = 'fieldmaps';
fieldmaps.name    = 'Input fieldmaps';
fieldmaps.help    = {['One fieldmap or 3 fieldmap gradient (dX dY dZ) files ' ...
                      'for optimizing BOLD sensitivity. Note that field ' ...
                      'derivatives must have units in T/m.']};
fieldmaps.ufilter = '.*';
fieldmaps.num     = [1 Inf];
% -------------------------------------------------------------------------
% template
% -------------------------------------------------------------------------
template         = cfg_files;
template.tag     = 'template';
template.name    = 'Input template';
template.help    = {'template for illustration (only placeholder at the moment)'};
template.ufilter = '.*';
template.num     = [1 Inf];
% -------------------------------------------------------------------------
% ROIs
% -------------------------------------------------------------------------
rois         = cfg_files;
rois.tag     = 'rois';
rois.name    = 'ROI files';
rois.help    = {'ROIs for which BOLD is optimized.'};
rois.ufilter = '.*';
rois.num     = [1 Inf];
% -------------------------------------------------------------------------
% Input Files
% -------------------------------------------------------------------------
inputfiles         = cfg_branch;
inputfiles.tag     = 'inputfiles';
inputfiles.name    = 'Input files';
inputfiles.val     = {fieldmaps template rois};
inputfiles.help    = {'Needed Input Files'};
% -------------------------------------------------------------------------
% menu main orientation
% -------------------------------------------------------------------------
main_orientation         = cfg_menu;
main_orientation.tag     = 'main_orientation';
main_orientation.name    = 'Choose the main orientation';
main_orientation.help    = {'Option to choose the main orientation'};
main_orientation.labels  = {'TRA' 'COR' 'SAG'};
main_orientation.values  = {'TRA' 'COR' 'SAG'};
main_orientation.val     = {'TRA'};
% -------------------------------------------------------------------------
% Field of View
% -------------------------------------------------------------------------
fov         = cfg_entry;
fov.tag     = 'fov';
fov.name    = 'Field of view';
fov.val     = {192};
fov.help    = {'Field Of View in Phase Encoding Direction in mm'};
fov.strtype = 'r';
fov.num     = [1 1];
% -------------------------------------------------------------------------
% Phase Resolution
% -------------------------------------------------------------------------
ph_res         = cfg_entry;
ph_res.tag     = 'ph_res';
ph_res.name    = 'Phase resolution';
ph_res.val     = {64};
ph_res.help    = {'Phase Resolution in #px'};
ph_res.strtype = 'r';
ph_res.num     = [1 1];
% -------------------------------------------------------------------------
% Phase Encoding Oversampling
% -------------------------------------------------------------------------
pe_ov         = cfg_entry;
pe_ov.tag     = 'pe_ov';
pe_ov.name    = 'Phase oversampling';
pe_ov.val     = {12};
pe_ov.help    = {'Oversampling Ratio in Phase Encoding Direction in %'};
pe_ov.strtype = 'r';
pe_ov.num     = [1 1];
% -------------------------------------------------------------------------
% Slice Thickness
% -------------------------------------------------------------------------
slicethickness         = cfg_entry;
slicethickness.tag     = 'slicethickness';
slicethickness.name    = 'Slice thickness';
slicethickness.val     = {3};
slicethickness.help    = {'Slice Thickness or (if known) the Full Width at ' ...
                          'Half Maximum (FWHM) of the slice excitation profile in mm'};
slicethickness.strtype = 'r';
slicethickness.num     = [1 1];
% -------------------------------------------------------------------------
% Echo Spacing
% -------------------------------------------------------------------------
echospacing         = cfg_entry;
echospacing.tag     = 'echospacing';
echospacing.name    = 'Echo spacing';
echospacing.val     = {0.5};
echospacing.help    = {'Echo Spacing in ms'};
echospacing.strtype = 'r';
echospacing.num     = [1 1];
% -------------------------------------------------------------------------
% Echo Time
% -------------------------------------------------------------------------
echotime         = cfg_entry;
echotime.tag     = 'echotime';
echotime.name    = 'Echo time';
echotime.val     = {30};
echotime.help    = {'Echo Time in ms'};
echotime.strtype = 'r';
echotime.num     = [1 1];
% -------------------------------------------------------------------------
% Voxel Size
% -------------------------------------------------------------------------
vox         = cfg_entry;
vox.tag     = 'vox';
vox.name    = 'Voxel size';
vox.val     = {[3 3 3]};
vox.help    = {'Voxel Size [read, phase, slice] in mm'}; % Note the order of input!
vox.strtype = 'r';
vox.num     = [1 3];
% -------------------------------------------------------------------------
% Acceleration Factor
% -------------------------------------------------------------------------
AccF         = cfg_entry;
AccF.tag     = 'AccF';
AccF.name    = 'Acceleration factor';
AccF.val     = {1};
AccF.help    = {'In-plane/GRAPPA Acceleration Factor'};
AccF.strtype = 'r';
AccF.num     = [1 1];
% -------------------------------------------------------------------------
% Partial Fourier Factor
% -------------------------------------------------------------------------
PF         = cfg_entry;
PF.tag     = 'PF';
PF.name    = 'Partial Fourier factor';
PF.val     = {1};
PF.help    = {'In-plane Partial Fourier Factor'};
PF.strtype = 'r';
PF.num     = [1 1];
% -------------------------------------------------------------------------
% Fixed Protocol Parameters
% -------------------------------------------------------------------------
fixedparameters         = cfg_branch;
fixedparameters.tag     = 'fixedparameters';
fixedparameters.name    = 'Fixed Protocol Parameters';
fixedparameters.val     = {main_orientation fov ph_res pe_ov slicethickness echospacing echotime vox AccF PF};
fixedparameters.help    = {'Fixed Protocol Parameters'};

% =========================================================================
% Set Simulation Parameters
% =========================================================================

% -------------------------------------------------------------------------
% Parameter shimz
% -------------------------------------------------------------------------
shimz         = cfg_entry;
shimz.tag     = 'shimz';
shimz.name    = 'shimz';
shimz.val     = {[-5 0 5 0.5]};
shimz.help    = {'Shim Gradient moment in z-direction in mT/m*ms [min ref max step-size]'};
shimz.strtype = 'r';
shimz.num     = [1 4];
% -------------------------------------------------------------------------
% Parameter tilt
% -------------------------------------------------------------------------
tilt         = cfg_entry;
tilt.tag     = 'tilt';
tilt.name    = 'tilt';
tilt.val     = {[-45 0 45 5]};
tilt.help    = {'tilt in degree [min ref max step-size]'};
tilt.strtype = 'r';
tilt.num     = [1 4];
% -------------------------------------------------------------------------
% Simulation Parameters
% -------------------------------------------------------------------------
simu         = cfg_branch;
simu.tag     = 'simu';
simu.name    = 'Simulation Parameters';
simu.val     = {shimz tilt};
simu.help    = {['Parameters to be simulated: all these parameters have a minimum, ' ...
                 'maximum and default value and a step size for the optimization procedure']};

% =========================================================================
% Other Settings
% =========================================================================

% -------------------------------------------------------------------------
% Reduce Field size
% -------------------------------------------------------------------------
rfs         = cfg_entry;
rfs.tag     = 'rfs';
rfs.name    = 'Reduce Field size';
rfs.val     = {0};
rfs.help    = {'0 = no (original size), 1 = yes (1/3)'};
rfs.strtype = 'r';
rfs.num     = [1 1];
% -------------------------------------------------------------------------
% Additional Inputs for R2s Option
% -------------------------------------------------------------------------
Global_3T         = cfg_entry;
Global_3T.tag     = 'Global_3T';
Global_3T.name    = '3T Global R2*';
Global_3T.val     = {1/45};                              
Global_3T.help    = {'Global R2* Value in 3T in [ms]^-1'};
Global_3T.num     = [1 1];
% -------------------------------------------------------------------------
% Additional Inputs for R2s Option
% -------------------------------------------------------------------------
Global_7T         = cfg_entry;
Global_7T.tag     = 'Global_7T';
Global_7T.name    = '7T Global R2*';
Global_7T.val     = {1/30};                              
Global_7T.help    = {'Global R2* Value in 7T in [ms]^-1'};
Global_7T.num     = [1 1];
% -------------------------------------------------------------------------
% Additional Inputs for R2s Option
% -------------------------------------------------------------------------
Voxel_wise         = cfg_files;
Voxel_wise.tag     = 'Voxel_wise';
Voxel_wise.name    = 'Voxel_wise R2* map';                             
Voxel_wise.help    = {'Input an R2* map in [ms]^-1'};
Voxel_wise.num     = [1 Inf];
% -------------------------------------------------------------------------
% Additional Inputs for R2s Option
% -------------------------------------------------------------------------
ROI_Averaged         = cfg_files;
ROI_Averaged.tag     = 'Voxel_wise';
ROI_Averaged.name    = 'Voxel_wise R2* for ROI-specific value';                            
ROI_Averaged.help    = {['Input an R2* map in [ms]^-1, an averaged value within' ...
                         'each ROI will be used for optimization']};
ROI_Averaged.num     = [1 1];
% -------------------------------------------------------------------------
% R2star Designation
% -------------------------------------------------------------------------
R2sOpt         = cfg_choice;
R2sOpt.tag     = 'R2sOpt';
R2sOpt.name    = 'R2star Option';
R2sOpt.values  = {Global_3T Global_7T Voxel_wise ROI_Averaged};
R2sOpt.val     = {Global_3T};
R2sOpt.help    = {'How to inoporate R2* in the optimization; The options are:' ...
                  '1. Global Value in 3T (1/45 ms^-1)' ...
                  '2. Global value in 7T (1/30 ms^-1)' ...
                  '3. Voxel-wise Map (ms^-1)' ...
                  '4. ROI-specific Averaged Value'
                  };
% -------------------------------------------------------------------------
% Other Settings
% -------------------------------------------------------------------------
other         = cfg_branch;
other.tag     = 'other';
other.name    = 'Other Settings';
other.val     = {rfs R2sOpt};
other.help    = {'Other Settings Used for Optimization'};

% =========================================================================
% Preprocessing
% =========================================================================
fmpoptbs         = cfg_exbranch;
fmpoptbs.tag     = 'FmpOptBS';
fmpoptbs.name    = 'BS Optimisation';
fmpoptbs.val     = {inputfiles fixedparameters simu other};
fmpoptbs.help    = {'This toolbox is currently only work in progress.'};
fmpoptbs.prog = @fmpoptbs_apply;
fmpoptbs.vout = @vout_fmpoptbs_apply;
%--------------------------------------------------------------------------
function opt = fmpoptbs_apply(job)

opt.results = epi_opt_param_TB(job.inputfiles.fieldmaps, job.inputfiles.rois, ...
                               job.inputfiles.template, ...
                               job.fixedparameters.main_orientation, ...
                               job.fixedparameters.fov, ...
                               job.fixedparameters.ph_res, ...
                               job.fixedparameters.pe_ov, ...
                               job.fixedparameters.slicethickness, ...
                               job.fixedparameters.echospacing, ...
                               job.fixedparameters.echotime, ...
                               job.fixedparameters.vox, ...
                               job.fixedparameters.AccF, job.fixedparameters.PF, ...
                               job.simu.tilt, job.simu.shimz, ...
                               job.other.rfs, job.other.R2sOpt, '_Opt');

% Not sure if this is necessary!
% opt.fmfiles = job.inputfiles.fieldmaps;

% =========================================================================
function dep = vout_fmpoptbs_apply(~)
% do something
dep = cfg_dep;

% end;
%--------------------------------------------------------------------------
