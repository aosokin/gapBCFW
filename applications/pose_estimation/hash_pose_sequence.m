function hash_code = hash_pose_sequence( y )

if size(y.mix_id{2},1)==3
    hash_code = [num2str((cat(1,y.mix_id{:}))') num2str(floor(y.joints(:))') ];
else
    hash_code = [num2str((cat(2,y.mix_id{:}))) num2str(floor(y.joints(:))') ];
end

end

