function result = epi_opt_param_TB(fieldmaps, rois, template, main_orientation, ...
                     fov, ph_res, pe_ov, delta_z, echo_spacing, TC, vx_epi, AF, ...
                     PF, tilt, shimz, rfs, R2sOpt, suffix)

% =========================================================================
% Copyright (C) 2015-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% =========================================================================

% ========================================================================= 
% This function iterates through the simulation parameter space to identify 
% the optimal BOLD (Blood Oxygen Level Dependent) sensitivity. 
% ========================================================================= 
% fieldmaps                       : Cell Array containing Field gradients 
%                                   in (read, phase, slice) directions in T/m          
% rois                            : Cell Array containing the ROIs
% template                        : Cell Array containing the Brain Mask       
% main_orientation                : Default Slice Orientation in EPI Acquistion                            
% fov                             : Field of View in the Phase Direction  (in mm)
% ph_res                          : Phase Resolution in number of pixels (Matrix size)
% pe_ov                           : Oversampling in Phase Encoding Direction in %
% delta_z                         : Slice Thickness or more precisely the
%                                   FWHM of the slice excitation profile
%                                   for a Guassian RF Pulse (in mm)
% echo_spacing                    : Echo spacing (in ms)
% TC                              : Central Echo Time (in ms)
% vx_epi                          : Voxel size (in mm)
% AF                              : In-plane Acceleration Factor
% PF                              : In-plane Partial Fourier
% tilt                            : Slice angulation (in degrees)        
%                                  (1x4) array [min ref max step-size]
% PP                              : Shim gradient moment in z-direction (in mT/m*ms)          
%                                  (1x4) array [min ref max step-size]  
% rfs                             : Reduced Field Size                     
%                                   0 = no (original size), 1 = yes (1/3) 
% R2sOpt                          : R2star Choice (in 1/ms)
% suffix                          : sufix for saved result
% =========================================================================
% Updated 28/09/2024
% by Shokoufeh Golshani

% -------------------------------------------------------------------------
% Define the output folder
% -------------------------------------------------------------------------
out_dir = fullfile(pwd, sprintf('results_BS%s', suffix));
if ~exist(out_dir,'file')
    mkdir(out_dir);
end

% -------------------------------------------------------------------------
% Reading Field Map files
% -------------------------------------------------------------------------
spm_progress_bar('Init', 10, 'preparing ...', 'steps');
spm_progress_bar('Set', 0);

for n = 1:length(fieldmaps)
    sprintf('Using fieldmap(s) = %s;',fieldmaps{n});
end

if (length(fieldmaps) == 3)
   fprintf('loading gradientmaps ...\n');

   rescale_gradient = 1.0;
   vol_fm_dX = spm_vol(fieldmaps{1});
   vol_fm_dY = spm_vol(fieldmaps{2});
   vol_fm_dZ = spm_vol(fieldmaps{3});

   fm_dX = resize(rescale_gradient.*spm_read_vols(vol_fm_dX), rfs);
   fm_dY = resize(rescale_gradient.*spm_read_vols(vol_fm_dY), rfs);
   fm_dZ = resize(rescale_gradient.*spm_read_vols(vol_fm_dZ), rfs);

elseif (length(fieldmaps) == 1)
   fprintf('loading fieldmaps and calculating gradientmaps ...\n');

   spm_progress_bar('Set', 1);
   vol_fm_dX = spm_vol(fieldmaps{1});

   [fm_dX, fm_dY, fm_dZ] = CalculateGradientmaps_TB(fieldmaps{1});
   fprintf('resizing gradientmaps ...\n');
   
   fm_dX = resize(fm_dX, rfs);
   fm_dY = resize(fm_dY, rfs);
   fm_dZ = resize(fm_dZ, rfs);

else
    fprintf('Error: invalid number of fieldmap files;\n');
end

% -------------------------------------------------------------------------
% Setting Protocol Parameters
% -------------------------------------------------------------------------
fprintf('setting user defined parameters ...\n');
spm_progress_bar('Set' ,7);

epi_param_fix.main_orientation = main_orientation;
epi_param_fix.echo_spacing = echo_spacing * 10^(-3);

epi_param_fix.fov      = fov * 10^(-3);
epi_param_fix.AccF     = AF;
epi_param_fix.ph_res   = ph_res;
epi_param_fix.delta_z  = delta_z * 10^(-3);
epi_param_fix.TC       = TC * 10^(-3);
epi_param_fix.vx_epi   = vx_epi * 10^(-3);

% Effective phase-encoding steps
epi_param_fix.pe_eff = ceil(epi_param_fix.ph_res * (1 + pe_ov/100));

% Fully-sampled case
epi_param_fix.TA_FS  = epi_param_fix.echo_spacing * epi_param_fix.pe_eff;

% PF-Acc case
epi_param_fix.pe_eff   = epi_param_fix.pe_eff * PF/epi_param_fix.AccF;

% Total acquisition time
epi_param_fix.TA       = epi_param_fix.echo_spacing * epi_param_fix.pe_eff; 

% -------------------------------------------------------------------------
% Setting Scanner-dependant Parameter
% -------------------------------------------------------------------------
fieldNames = fieldnames(R2sOpt);
R2sfield = fieldNames{1};

R2s = R2sOpt.(R2sfield);

if isnumeric(R2s)
    scanner_param.R2s = R2s;
    fprintf('loading R2s value: %0.2f (ms^-1) \n', R2s);
else
    vol_R2s = spm_vol(char(R2s));
    scanner_param.R2s = spm_read_vols(vol_R2s);
    fprintf('loading R2s map: %s\n', char(R2s));
end

scanner_param.R2s = scanner_param.R2s .* 10^3;

% -------------------------------------------------------------------------
% Reading ROIs
% -------------------------------------------------------------------------
fprintf('reading ROIs ...\n');
spm_progress_bar('Set', 8);

for n = 1:length(rois)
    sprintf('Reading ROI: %s;', rois{n});
    vol_MyRoi = spm_vol(rois{n});
    ROI_slct(:,:,:,n) = resize(spm_read_vols(vol_MyRoi), rfs);

    % ROI_averaged Susceptiblity Gradient
    GSSroi = fm_dZ(squeeze(ROI_slct(:,:,:,n))>0);
    GSSroi = mean(GSSroi);
    GSSroi_TE = GSSroi*epi_param_fix.TC*1e6;

    sprintf(['Mean slice gradient moment in the roi: %0.2f (mT/m*ms). ' ...
             'You can adjusct shimz range accordingly\n'], GSSroi_TE);

    % ROI_averaged R2s
    if strcmp(R2sfield, 'ROI_Averaged')
        R2sroi = scanner_param.R2s.*squeeze(ROI_slct(:,:,:,n));
        R2sroi = mean(nonzeros(R2sroi));
        scanner_param.R2s = R2sroi;
        fprintf('ROI averaged R2s value: %0.2f (s^-1) \n', R2sroi);
    end
end

% -------------------------------------------------------------------------
% Reading the template brain mask
% -------------------------------------------------------------------------
fprintf('reading the template brain mask ...\n');
spm_progress_bar('Set', 9);

for n = 1:length(template)
    vol_tmpl_msk = spm_vol(template{n});
    Brainmask_tmpl(:,:,:,n) = resize(spm_read_vols(vol_tmpl_msk), rfs);

end

% -------------------------------------------------------------------------
% Setting Simulation Parameters
% -------------------------------------------------------------------------
PP_range = shimz(1):shimz(4):shimz(3);
PP_ref = shimz(2);

tilt_range = tilt(1):tilt(4):tilt(3);
tilt_ref = tilt(2); 

% =========================================================================
% Phase Encoding direction (based on the prewinder gradient moment)
% (-1 = Positive Prewinder, +1 = Negative Prewinder)
% =========================================================================
PE_range = [-1 1];

% -------------------------------------------------------------------------
% Setting Field derivatives
% -------------------------------------------------------------------------
FG.DX        = fm_dX;
FG.DY        = fm_dY;
FG.DZ        = fm_dZ;
FG.direction = vol_fm_dX.mat(1:3,1:3);

% -------------------------------------------------------------------------
% BS for Baseline Protocol (no tilt, no compensation, Positive Prewinder (AP))
% -------------------------------------------------------------------------
epi_param_opt.GP = [0 0 PP_ref]*10^-6;              % Compensation gradient moment (T/m*s)
epi_param_opt.tilt = tilt_ref;                      % Tilt of slice (deg, T>C, Siemens convention)
epi_param_opt.PE_dir = PE_range(1);                 % PE direction of EPI 

BS_baseline = CalculateBS_TB(FG, epi_param_opt, epi_param_fix, scanner_param);

% -------------------------------------------------------------------------
% Exploring the Parameter Space
% -------------------------------------------------------------------------
fprintf('exploring parameter space ... \n');

result.BS_matrix = zeros(2, size(tilt_range, 2), size(PP_range, 2), length(rois));
size_parameterspace = 2*size(tilt_range, 2)*size(PP_range, 2);

spm_progress_bar('Clear');
spm_progress_bar('Init', size_parameterspace, 'exploring parameter space', 'settings completed');

ct = 0;
ct0 = 0;

for PE_val = 1:length(PE_range)
    
    for tilt_val = 1:length(tilt_range)
        
        for PP_val = 1:length(PP_range)
            
            spm_progress_bar('Set', ct);
            ct = ct + 1;

            PE_all(ct) = PE_range(PE_val);
            tilt_all(ct) = tilt_range(tilt_val);
            PP_all(ct) = PP_range(PP_val);

            if (PE_range(PE_val) == -1) && (tilt_range(tilt_val) == 0) && (PP_range(PP_val) == 0)
                ct0 = ct;     
            end

            epi_param_opt.GP = [0 0 PP_range(PP_val)]*10^-6;
            epi_param_opt.tilt = tilt_range(tilt_val);
            epi_param_opt.PE_dir = PE_range(PE_val);
      
            BS = CalculateBS_TB(FG, epi_param_opt, epi_param_fix, scanner_param);
      
            % -------------------------------------------------------------
            % Excluding out of range values
            % -------------------------------------------------------------
            BS_gain = (BS./(BS_baseline+eps)-1)*100;
            BSGainMask = (BS_gain > -100) & (BS_gain < 200);
            BrainAndGainMask = (Brainmask_tmpl > 0.99).*BSGainMask;

            % -------------------------------------------------------------
            % Evaluating ROIS
            % -------------------------------------------------------------
            for n = 1:length(rois)
                Ind_ToOpt = find((BrainAndGainMask.*squeeze(ROI_slct(:,:,:,n))) > 0);
                Roi_Sel_val = BS(Ind_ToOpt);
                Roi_mean(n,ct) = mean(Roi_Sel_val(:));
                Roi_std(n,ct)  = std(Roi_Sel_val(:));
                result.BS_matrix(PE_val, tilt_val, PP_val, n) = Roi_mean(n, ct);
            end

        end
    end
end

spm_progress_bar('Clear');

fprintf('------------------------------------------------------------------\n');
fprintf('BS Optimization done for\n');

for n = 1:length(rois)
    display(sprintf('ROI Nr. %2d: %s;',n,rois{n}));
end

fprintf(' \n');
fprintf('Optimal parameters:\n');

for n = 1:length(rois)
   [v, I] = max(squeeze(Roi_mean(n,:)));
   fprintf('ROI Nr.: %2d; BS: %6.3f; BS-REF: %6.3f; BS-gain: %6.3f; PE: %1d; PP: %4.1f; tilt: %4d;\n', ...
            n, v, Roi_mean(n,ct0), ((v/Roi_mean(n,ct0))-1)*100, PE_all(I), PP_all(I), tilt_all(I));
   result.results(n, 1) = I;                       result.results(n, 2) = v; 
   result.results(n, 3) =  v/Roi_mean(n,ct0);      result.results(n, 4) = PE_all(I); 
   result.results(n, 5) = PP_all(I);               result.results(n, 6) = tilt_all(I);
   result.results(n, 7) = Roi_mean(n, ct0);        result.results(n, 8) = Roi_std(n,I);
end

fprintf(' \n');
fprintf('BS = Optimized relative BS in ROI \n');
fprintf('BS-REF = Baseline relative BS in ROI \n');
fprintf('BS-gain = BS in optimal protocol compared to default protocol\n');
fprintf('PE = phase encoding direction (-1 = AP, +1 = PA)\n');
fprintf('PP = Shim gradient moment in z-direction (mT/m*ms)\n');
fprintf('tilt = tilt of slice (Positive = towards foot, Negative = towards head)\n');
fprintf('------------------------------------------------------------------\n');

% -------------------------------------------------------------------------
% Saving the result
% -------------------------------------------------------------------------
matname = fullfile(out_dir, 'BSOpt.mat');
save(matname, 'result', 'BS', 'BS_baseline');

% -------------------------------------------------------------------------
% Display the result
% -------------------------------------------------------------------------
DisplayBS_TB(result, tilt_range, PP_range, rois)

result = 1;

end
