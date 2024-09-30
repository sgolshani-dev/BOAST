
toolboxFileLoc = [fullfile(spm('Dir'), 'toolbox','FmpOptBS') filesep];
jobfile = {[fullfile(toolboxFileLoc, 'BS_Optimisation_job.m')]};

nrun = X;                                  % enter the number of runs here
jobs = repmat(jobfile, 1, nrun);
inputs = cell(0, nrun);
for crun = 1:nrun
end
spm('defaults', 'FMRI');
spm_jobman('run', jobs, inputs{:});