% function saveCytosimOutput(simPath, nModes, savefolder)
% Function to read output from JFilament and save filament data, filament lengths,
% and bending mode amplitudes in savefolder as csv files. Make sure to run this
% after running cytosimDataExtractor.m, found in this same folder.
%
% Parameters
% ----------
% simPath : str
%     location of simulation output data after reformating using extractor.m
% nModes : int
%     highest mode amplitude to calculate
% px2um : px2um
%     microns per pixel from microscopy images
% savefolder : str
%     string of absolute path to save output csv files.
%
% Returns
% -------
% savefolder/filamentLengths.csv
% savefolder/modeCoeffs.csv
%
% Created by Daniel Seara 10/20/2018
function saveCytosimOutput(simPath, nModes, savefolder)
    pathParts = split(simPath, filesep);
    disp(['Reading data from ', pathParts{end}, '...'])
    simData = dlmread(fullfile(simPath, [pathParts{end}, '.dat']), ',');
    disp('done')

    fids = unique(simData(:, 3));
    frames = unique(simData(:, 4));

    isOutside = zeros(size(fids));
    for ii = 1:numel(fids)
        data = simData(simData(:,3) == fids(ii), [1,2]);
        isOutside(ii) = checkFilPosition(data(:,1), data(:,2), 10, 10, 0.95);
    end

    insideFids = fids(~isOutside);

    simDataInside = simData(ismember(simData(:,3), insideFids), :);

    tic
    disp('Calculating modes...')
    [aaAll, avgL] = getModes(simDataInside, nModes, 'free');
    disp('done')
    toc

    aaAll = [aaAll, repmat(frames, numel(insideFids), 1), repelem(insideFids, numel(frames))];
    avgL = [avgL, insideFids];

    dlmwrite(fullfile(savefolder, 'modeCoeffs.csv'), aaAll, ',')
    dlmwrite(fullfile(savefolder, 'filamentLengths.csv'), avgL, ',')
