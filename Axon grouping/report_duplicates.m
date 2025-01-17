function duplicates = report_duplicates(groups)
    memberDict = containers.Map('KeyType', 'double', 'ValueType', 'any');
    
    for i = 1:length(groups)
        for j = 1:length(groups{i})
            member = groups{i}(j);
            if memberDict.isKey(member)
                memberDict(member) = [memberDict(member), i];
            else
                memberDict(member) = [i];
            end
        end
    end

    duplicates = {};
    for key = keys(memberDict)
        loc = memberDict(key{1});
        if length(loc) > 1
            fprintf('Member %d appears in groups: %s\n', key{1}, mat2str(loc));
            duplicates{end+1} = struct('Member', key{1}, 'Groups', loc);
        end
    end
    
    if isempty(duplicates)
        fprintf('No duplicates found.\n');
    end
end
