function build_distance_transform()

src_path = fullfile('helpers', 'mex');
src_files = {'distance_transform.cpp'};
output_path = 'helpers';

mexcmd = ['mex -O ', '-outdir ', output_path];
for i_file = 1 : length(src_files)
    eval([mexcmd,  ' ', fullfile(src_path, src_files{i_file})]);
end
