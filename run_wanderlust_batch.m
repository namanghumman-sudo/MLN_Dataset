try
    cd(fileparts(mfilename('fullpath')));
    addpath(genpath(pwd));
    run('run_wanderlust_export.m');
catch ME
    disp(getReport(ME, 'extended'));
    exit(1);
end
exit(0);
