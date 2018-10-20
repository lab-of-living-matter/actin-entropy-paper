% function saveAxonemeOutput(axonemeFile, nModes, px2um, dt, savefolder)
%
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
% dt : float
%     time between frames. Used to subtract linear increase of tangent angle
%     measured as axoneme rotates continuously
% savefolder : str
%     string of absolute path to save output csv files.
%
% Returns
% -------
% savefolder/filamentLength.csv
% savefolder/modeCoeffs.csv
%
% Created by Daniel Seara 10/20/2018
function saveAxonemeOutput(axonemeFile, nModes, px2um, dt, savefolder)
    axoneme = dlmread(axonemeFile, '\t', 1, 0);
    aa = getAxonemeModes(axoneme, nModes, px2um, dt);
    L = px2um .* size(aa, 1);
    dlmwrite(fullfile(savefolder, 'filamentLength.csv'), L, ',')
    dlmwrite(fullfile(savefolder, 'modeCoeffs.csv'), aa, ',')

