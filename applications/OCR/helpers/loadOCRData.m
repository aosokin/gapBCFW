function [patterns_train, labels_train, patterns_test, labels_test] = loadOCRData(data_name, data_path)

% load matlab-compatible OCR data, if it doesn't exist, create the .mat file
fname = fullfile(data_path, 'ocr.mat');
if (~exist(fname, 'file'))
    fprintf('Creating ocr.mat for the first time...\n')
    convertOCR(data_path);
else
    fprintf('Loading already built ocr.mat file...\n')
end
temp = load(fname);
dataset = temp.dataset;

% split data into training and testing chuncks and construct
% the x_i data (patterns) and y_i labels cell arrays. 
patterns_train = {}; % for training
labels_train = {};
patterns_test = {}; % for testing
labels_test = {};
idx_train = 1;
idx_test = 1;
for i=1:length(dataset)
    
    % number of states an individual y_i can assume
    num_states = 26;
    [xi,yi] = constructExampleOCR(dataset{i}, num_states);
    
    if ( (isequal(data_name, 'ocr2')  && dataset{i}.fold ~= 0) ...
        || (isequal(data_name, 'ocr') && dataset{i}.fold == 0))
        patterns_train{idx_train} = xi;
        labels_train{idx_train} = yi;
        idx_train = idx_train+1;
    else
        patterns_test{idx_test} = xi;
        labels_test{idx_test} = yi;
        idx_test = idx_test+1;
    end
end
fprintf('successfully loaded %s dataset.\n', data_name)

end % loadOCRData


function [pattern,label] = constructExampleOCR(datapoint, num_states)

pattern = [];
label = datapoint.word; % starts with 0
pattern.data = [];
pattern.num_states = num_states;
for p_idx=1:numel(datapoint.pixels)
    t = double(datapoint.pixels{p_idx}(:));
    pattern.data = [pattern.data [t; 1.0]];
end

end % constructExampleOCR


function convertOCR(data_path)

fname_data = fullfile(data_path, 'letter.data');
fname_names = fullfile(data_path, 'letter.names');
if ~exist(fname_data, 'file') || ~exist(fname_names, 'file')
    fprintf('Files\n%s\nor\n%s\nnot found. Trying to download...\n', fname_data, fname_names);
    url_letter_names = 'http://www.seas.upenn.edu/~taskar/ocr/letter.names';
    url_letter_data = 'http://www.seas.upenn.edu/~taskar/ocr/letter.data.gz';
    urlwrite(url_letter_names, fname_names);
    gunzip(url_letter_data, data_path);
    if ~exist(fname_data, 'file') || ~exist(fname_names, 'file')
        fprintf('Files\n%s\nor\n%s\n could not be downloaded. Check the internet connection or download them manually from\n%sand%s\n', fname_data, fname_names, url_letter_data, url_letter_names);
    else
        fprintf('Download complete. Processing...\n');
    end
end
fid = fopen(fname_data);
c = textscan(fid,'%d%s%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d', 'Delimiter', '\t');
fclose(fid);

dataset = {};
first = 1;
P = cell2mat(c(7:end));
for i=1:length(c{1})
    if (first)
        d = [];
        d.fold = c{6}(i);
        d.pixels = {};
        d.word = '';
        first = 0;
    end

    % append
    d.word = sprintf('%s%s', d.word, c{2}{i});
    d.pixels = [d.pixels reshape(P(i, :), [8 16])'];

    if (c{3}(i) == -1)
        dataset = [dataset; d];
        first = 1;
    end
end

for i=1:length(dataset)
    d = dataset{i};
    temp = length(d.word);
    for j=1:length(d.word)
        temp(j) = d.word(j) - 'a';
    end
    d.word = temp;
    dataset{i} = d;
end

fname = fullfile(data_path, 'ocr.mat');
save(fname, 'dataset');

end % convertOCR
