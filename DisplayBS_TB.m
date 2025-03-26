function DisplayBS_TB(result, tilt_range, zshim_range, rois)

nROIs = length(rois); 

BS_matrix = result.BS_matrix;
BS_baseline = result.results(1:nROIs, 7);

%==========================================================================
% We assume only 15% drop in wellshimmed areas is acceptable (BS = 1)
%==========================================================================
minBS = 0.85;             
maxBS = 1.3;

% Determine optimal figure size based on number of ROIs
figureSize = max(800, 100 * nROIs);          % Increase figure size dynamically
numRows = ceil(nROIs / 5) * 2;               % Adjust rows for large nROIs
numCols = min(nROIs, 5);                     % Keep max 5 columns per row

figure(1); 
set(gcf, 'Color', [1 1 1], 'Position', [100, 100, figureSize, figureSize/2]);

for ROIset = 1:nROIs
    % Display first BS matrix for each ROI
    subplot(numRows, numCols, ROIset);
    imagesc(zshim_range, tilt_range, squeeze(BS_matrix(1, :, :, ROIset)), [minBS maxBS]);
    title(spm_file(rois{ROIset}, 'filename'), 'Interpreter', 'none');
    xlabel('Z-Shim');
    ylabel('Tilt');
    colorbar;
    
    % Display second BS matrix for each ROI
    subplot(numRows, numCols, ROIset + numCols);
    imagesc(zshim_range, tilt_range, squeeze(BS_matrix(2, :, :, ROIset)), [minBS maxBS]);
    title(spm_file(rois{ROIset}, 'filename'), 'Interpreter', 'none');
    xlabel('Z-Shim');
    ylabel('Tilt');
    colorbar;
end
sgtitle('BOLD Sensitivity Matrices for All ROIs'); % Super title for clarity

% ------------------------ BOLD Sensitivity Stats ------------------------
figure(5); 
set(gcf, 'Color', [1 1 1], 'Position', [100, 100, 800, 500]);

subplot(2, 3, 1)
bar(result.results(:, 2));
hold on; 
bar(BS_baseline);
hold off;
title('BOLD Sensitivity of Optimum');

subplot(2, 3, [2 3])
bar(result.results(:, 3)); 
title('Relative BOLD Sensitivity Gain of Optimum')

subplot(2, 3, 4)
bar(result.results(:, 4));
title('PE Polarity')

subplot(2, 3, 5)
bar(result.results(:, 5));
title('Z shim')

subplot(2, 3, 6)
bar(result.results(:, 6));
title('Tilt')
