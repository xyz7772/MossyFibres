function [is_independent, detected_indices] = check_MF_independence(b, idx, group_ids)
    % b: Input vector, containing the elements to be checked
    % group_ids: A cell array containing multiple groups, each group is an array
    % is_independent: Returned boolean value, indicating whether the elements are fully independent
    % detected_indices: Indices of the elements that are not independent

    % Extract the first four elements from b
    elements = b(idx)';
    detected_indices = [];

    has_common_mf = false;

    for i = 1:length(group_ids)
        current_group = group_ids{i};
        current_group_with_index = [current_group, i];

        % Check if any element in b is present in the current group
        if sum(ismember(elements, current_group_with_index)) >= 2
            has_common_mf = true;
            detected_indices = find(ismember(elements, current_group_with_index) == 1);
            fprintf('MFB %s found in %s, belong to the same MF\n', mat2str(b(detected_indices)),  mat2str(current_group_with_index));
            break;
        end
    end

    if has_common_mf
        is_independent = false; % Not fully independent
    else
        is_independent = true; % Fully independent
        fprintf('MFBs are independent\n')
    end
end
