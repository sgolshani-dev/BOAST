function scanner_param = SetDefaultScannerParam
% This is no longer needed!
% =========================================================================
% This function sets the scanner parameters. 
% These values can be modified as needed to suit specific requirements.
% =========================================================================
% name                            : Scanner name
% B0                              : Field strength
% T2s                             : T2* value (in ms) - An averaged value over 
%                                   the entire brain or a voxel-wise map.
%                                   Although not directly a scanner parameter, 
%                                   it depends on the field strength!
% =========================================================================

% Updated 23/09/2024
% by Shokoufeh Golshani

% This folder is not needed!
scanner_param.name = 'Trio';
scanner_param.B0 = 3;

% well this is not directly a scanner parameter, however it depends on the
% field strength and therefore scanner
scanner_param.T2s = 45; 	    % T2* at 3 T (Wansapura et al., JMRI 1999)


end
