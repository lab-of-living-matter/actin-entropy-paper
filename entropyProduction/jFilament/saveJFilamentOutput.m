% function saveJFilamentOutput(filaments, nModes, px2um, savefolder)
% Function to read output from JFilament and save filament data, filament lengths,
% and bending mode amplitudes in savefolder as csv files.
%
% Parameters
% ----------
% filaments : matlab table
%     result of importing txt file that JFilament saves by default. Replace empty
%     elements with NaNs
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
% savefolder/filamentData.csv
%
% Created by Daniel Seara 10/20/2018
function saveJFilamentOutput(filaments, nModes, px2um, savefolder)
    filaments_parsed = loadJFilamentData(filaments);
    [aa, aaEach, avgL] = getJFilamentModes(filaments_parsed, nModes, px2um);
    filaments_parsed_onematrix = vertcat(filaments_parsed{:});
    filaments_parsed_onematrix(:,end) = repelem(1:numel(filaments_parsed), cellfun(@(x) (size(x,1)), filaments_parsed))';
    aa = [aa, repelem(1:numel(aaEach), cellfun(@(x) (size(x,1)), aaEach))'];

    dlmwrite(fullfile(savefolder, 'filamentLengths.csv'), avgL, ',')
    dlmwrite(fullfile(savefolder, 'modeCoeffs.csv'), aa, ',')
    dlmwrite(fullfile(savefolder, 'filamentData.csv'), filaments_parsed_onematrix, ',')
