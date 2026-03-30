% call_run_tspace.m
here = fileparts(mfilename('fullpath'));
cd(here);

logFile = fullfile(here, 'matlab_tspace_log.txt');
diary(logFile); diary on;

disp(['PWD: ', pwd]);
disp('Calling run_tspace...');

try
    run_tspace;                 % calls run_tspace.m in this folder
    disp('run_tspace finished OK');
catch ME
    disp('*** MATLAB ERROR ***');
    disp(getReport(ME,'extended'));
    diary off;
    exit(1);
end

diary off;
exit(0);