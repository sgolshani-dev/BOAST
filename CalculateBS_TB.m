function BS = CalculateBS_TB(FG, epi_param_opt, epi_param_fix, scanner_param) 

% ========================================================================
% This function calculates BOLD sensitivity using field map gradients and 
% a defined set of parameters.
%% calculate BOLD sensitivity from gradient fieldmap
%% based on calc_BS_fm_atlas (NW)    % Not sure what this atlas is
% Copyright (C) 2014-2018 Steffen Volz
% Wellcome Trust Centre for Neuroimaging, London
% and Max Planck Institute for Human Cognitive and Brain Sciences, Leipzig 
% ========================================================================

% Updated 28/09/2024
% By Shokoufeh Golshani

% =========================================================================
% Unpack Input Variables
% =========================================================================
gam = 2*pi*42.58e6;                     % gyromagnetic ratio for protons in Hz/T

default_epi_params = SetDefaultEPIParam;
try fov = epi_param_fix.fov;     catch,   fov = default_epi_params.fov;         end
try AcF = epi_param_fix.AccF;    catch,   AcF = default_epi_params.AccF;        end
try PF = epi_param_fix.PF;       catch,   PF = default_epi_params.PF;           end
try TC = epi_param_fix.TC;       catch,   TC = default_epi_params.echotime;     end
try TA = epi_param_fix.TA;       catch,   TA = default_epi_params.TA;           end
try TA_FS = epi_param_fix.TA_FS; catch,   TA_FS = default_epi_params.TA_FS;     end
try vx_epi = epi_param_fix.vx_epi; catch, vx_epi = default_epi_params.vox;      end
try delta_z = epi_param_fix.delta_z; ...
catch, delta_z = default_epi_params.delta_z;   end
try echo_spacing = epi_param_fix.echo_spacing; ...      
catch, echo_spacing = default_epi_params.echo_spacing;     end
try main_orientation = epi_param_fix.main_orientation;     
catch,   main_orientation = default_epi_params.main_orientation;      end

% Field gradients
fm_dX = FG.DX;
fm_dY = FG.DY;
fm_dZ = FG.DZ;
direction   = FG.direction;

% =========================================================================
% Compensation Gradient Moments
% =========================================================================
GPrep_RO = epi_param_opt.GP(1);
GPrep_PE = epi_param_opt.GP(2);
GPrep_SS = epi_param_opt.GP(3);

% =========================================================================
% Rotation matrix for converting filed map gradients from XYZ to RPS
% =========================================================================
Angle = epi_param_opt.tilt/180*pi;

% RO: RL, PE: PA, SL: SI;  Rotation about the RO axis (Positive towards foot)
if strcmp(main_orientation,'TRA') == 1
    Rrot = [1     0           0;
            0  cos(Angle) -sin(Angle);
            0  sin(Angle)  cos(Angle)];

% Transformation from RAS+ (Nifti format) to LAI+ (Scanner Coordinte)
% this can be translated as a 180-degree rotation around the second axis
Rtrans = [-1  0  0;
           0  1  0;
           0  0 -1]*direction;

% RO: PA, PE: IS, SL: LR;  Rotation about the RO axis    
elseif strcmp(main_orientation,'SAG') == 1
    Rrot = [0  -sin(Angle)  cos(Angle);
            0   cos(Angle)  sin(Angle);
            -1     0           0];

% RO: LR, PE: SI, SL: PA;  Rotation about the RO axis   
elseif strcmp(main_orientation,'COR') == 1
    Rrot = [-1     0           0;
            0  -sin(Angle)  cos(Angle);
            0   cos(Angle)  sin(Angle)];
end

Rtotal = Rrot*Rtrans;

fGR = fm_dX*Rtotal(1,1) + fm_dY*Rtotal(1,2) + fm_dZ*Rtotal(1,3);
fGP = fm_dX*Rtotal(2,1) + fm_dY*Rtotal(2,2) + fm_dZ*Rtotal(2,3);
fGS = fm_dX*Rtotal(3,1) + fm_dY*Rtotal(3,2) + fm_dZ*Rtotal(3,3);

% =========================================================================
% Calculate the Q value which determines distortion and echo shift
% =========================================================================
if epi_param_opt.PE_dir == 1                  % A negative prephasing and positive blips          
    Q = 1 + (gam * echo_spacing * (fov/AcF) * fGP);
else                                          % A positive prephasing and negative blips 
    Q = 1 - (gam * echo_spacing * (fov/AcF) * fGP);
end

% =========================================================================
% Actual (local) Echo time
% =========================================================================
% if no shimming in the PE direction is used, this will be zero!
TE_shift = gam * GPrep_PE * (fov/AcF) * echo_spacing;

TE = (TC + TE_shift)./Q;
dTE = (TE - TC); 

% =========================================================================
% Sudden dropout when echo is shifted out of Acquisition window in the PE 
% and RO directions
% =========================================================================
shift_mask = ones(size(Q));

if PF ~= 1
    shift_mask((dTE) < -(TA - TA_FS/2/AcF)) = 0;
    shift_mask((dTE) > (TA_FS/2/AcF)) = 0;
else
    shift_mask(abs(dTE) > (TA/2)) = 0;
end

shift_mask = shift_mask.*(abs(fGR + GPrep_RO./TE) < (pi/gam./TE/vx_epi(1)));

% =========================================================================
% Contribution of the through-plane field gradient --- Gaussian RF pulse
% =========================================================================
I = exp((-gam^2*delta_z^2/16/log(2)).*((GPrep_SS + fGS.*TE).^2));


% =========================================================================
% BOLD Sensitivity
% =========================================================================
BS = (I./Q.^2).*exp(-(TC.*scanner_param.R2s).*((1./Q)-1));

% =========================================================================
% Correct maps for shifts of data out of acquisition window
% =========================================================================
BS = BS.*shift_mask;

end

