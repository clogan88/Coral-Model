function [multiThread, queueMax] = parallelSetup(n)
    % Cases: n undefined - use existing pool
    %        n = existing pool size - use exising pool
    %        n = 0 - leave existing pool alone (if any). Don't use it.
    %        n != existing pool size - delete any old pool and create one
    if nargin == 0
        % User wants to use whatever exists.
        pool = gcp('nocreate'); % Maximum number of threads for PDF creation.
        if isempty(pool)
            threads = 0;
            multiThread = false;
        else
            threads = pool.NumWorkers;
            multiThread = true;
        end
        fprintf('Using existing pool size with %d threads.\n', threads);
        poolReady = true;
    elseif n == 0
        % User wants no parallel threads.
        threads = 0;
        multiThread = false;
        poolReady = true;
    else
        % User specified a number of threads.  This will fail if more are
        % asked for than the value in Parallel Preferences.
        pool = gcp('nocreate');
        if isempty(pool)
            threads = n;
            multiThread = true;
            poolReady = false;
        elseif n == pool.NumWorkers
            threads = n;
            multiThread = true;
            poolReady = true;
        else
            threads = n;
            multiThread = true;
            poolReady = false;
            delete(gcp('nocreate'));
        end 
    end
    % Start a pool if needed, and notify in most cases.
    if multiThread && ~poolReady
        pool = parpool(threads);
        fprintf('Running with up to %d print threads.\n', pool.NumWorkers);
    elseif nargin == 0 && ~multiThread
        fprintf('Multithreaded PDF generation is off.  Start a parallel pool to enable it.');
    end
    queueMax = threads;
end