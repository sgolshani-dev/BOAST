function epi_param = SetDefaultEPIParam

% =========================================================================
% This function sets the fixed parameters for the EPI protocol. 
% These values can be modified as needed to suit specific requirements 
% or preferences.
% =========================================================================
% main_orientation                : Slice oriantation
%                                   'TRA' : transverse 
%                                   'CRO' : coronal
%                                   'SAG' : sagittal
% fov                             : Field of view (in mm)
% ph_res                          : Basic resolution in PE direction (Matrix size)
% pe_ov                           : Oversampling Ratio in Phase Encoding 
%                                   Direction in %
% slicethickness                  : Full width at half-maximum (FWHM) 
%                                   of the slice profile (in mm)
% echo_spacing                    : Echo spacing (in ms)
% echotime                        : Effective (central) echo time (in ms)
% vox                             : Voxel size (in mm) 
%                                   (1x3) array (read, phase, slice) direction
% AccF                            : In-plane Acceleration Factor
% PF                              : Partial Fourier Coefficient
% =========================================================================

% Updated 23/09/2024
% by Shokoufeh Golshani

epi_param.main_orientation = 'TRA';
epi_param.fov              = 192;
epi_param.ph_res           = 64;
epi_param.pe_ov            = 12;
% Note here 2 mm is used as the FWHM
epi_param.slicethickness   = 2; % Siemens pulse approximates Gaussian with 2 mm FWHM
epi_param.echo_spacing     = 0.5; 
epi_param.echotime         = 30; 
epi_param.vox              = [3 3 3];
epi_param.AccF             = 1;
epi_param.PF               = 1;


epi_param.fov              = epi_param.fov * 10^(-3);
epi_param.echotime         = epi_param.echotime * 10^(-3);
epi_param.echo_spacing     = epi_param.echo_spacing * 10^(-3);
epi_param.slicethickness   = epi_param.slicethickness * 10^(-3);
epi_param.vox              = epi_param.vox * 10^(-3);

% Effective phase-encoding steps
epi_param.pe_eff = ceil(epi_param.ph_res * (1 + epi_param.pe_ov/100));

% Fully-sampled case
epi_param.TA_FS  = epi_param.echo_spacing * epi_param.pe_eff;

% PF-Acc case
epi_param.pe_eff   = epi_param.pe_eff * epi_param.PF/epi_param.AccF;

% Total acquisition time
epi_param.TA       = epi_param.echo_spacing * epi_param.pe_eff;   


end
