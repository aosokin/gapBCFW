function hash_code = hash_segmentation(label_struct)
% Create a hash_code for horseSeg
% hash_code is simply string of '0' and '1' for the vector of labels
    label_array = label_struct.labels;
    hash_code = repmat('0', 1,length(label_array));
    hash_code(label_array==1) = '1';
end

