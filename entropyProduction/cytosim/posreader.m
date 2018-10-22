% Initialize variables.
function posfile = posreader(file)

    filename = file; % full path name
    delimiter = ' ';
    formatSpec = '%s%s%s%s%s%[^\n\r]';

    % Open the text file.
    fileID = fopen(filename,'r');

    % Read columns of data according to the format.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);

    % Close the text file.
    fclose(fileID);

    % Create output variable
    posfile = [dataArray{1:end-1}];

end
