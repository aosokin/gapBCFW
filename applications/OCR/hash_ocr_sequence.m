function hash_code = hash_ocr_sequence(label_array)
% Create a hash_code for the ocr sequence
% hash_code is simply string of characters:
% 0 = 'a', 1 = 'b', etc.
    alphabet = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z'};
    hash_code = repmat('a', 1,length(label_array));
    for i = 1:length(label_array)
        hash_code(i) = alphabet{label_array(i)+1};
    end
end

