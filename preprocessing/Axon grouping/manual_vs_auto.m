function [TP, TN, FP, FN] = manual_vs_auto(merged_manual, merged_ldr, merged_cc)

    TP = 0;
    TN = 0;
    FP = 0;
    FN = 0;

    ldr_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
    maunal_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    for j = 1:length(merged_ldr)
        primary = merged_ldr{j}(1);
        j
        ldr_map(primary) = merged_ldr{j};
    end
    
    for i = 1:length(merged_manual)
        primary = merged_manual{i}(1);
        i
        maunal_map(primary) = merged_manual{i};
    end

    for i = 1:length(merged_cc)
        cc_group = merged_cc{i};
        primary = cc_group(1);
        
        inLdr = ldr_map.isKey(primary);
        inGroups = maunal_map.isKey(primary);
        
        for member = cc_group
            ldr_presence = inLdr && any(member == ldr_map(primary));
            manual_presence = inGroups && any(member == maunal_map(primary));
            if ldr_presence && manual_presence
                TP = TP + 1;
            elseif ~ldr_presence && ~manual_presence
                TN = TN + 1;
            elseif ldr_presence && ~manual_presence
                FP = FP + 1; 
            elseif ~ldr_presence && manual_presence
                FN = FN + 1; 
            end
        end
    end
end
