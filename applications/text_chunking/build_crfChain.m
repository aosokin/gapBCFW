function build_crfChain()

src_path = fullfile('helpers', 'mex');
src_files = {'crfChain_makeLogNodePotentialsC.c', 'crfChain_maxMarginUpdateGradient.c'};
output_path = 'helpers';

mexcmd = ['mex -O ', '-outdir ', output_path];
for i_file = 1 : length(src_files)
    eval([mexcmd,  ' ', fullfile(src_path, src_files{i_file})]);
end
