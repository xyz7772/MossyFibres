function category = getMice(fileName, Jeremy_names, Bernie_names, Bill_names, Nigel_names)
    if ismember(fileName, Jeremy_names)
        category = 'Jeremy';
    elseif ismember(fileName, Bernie_names)
        category = 'Bernie';
    elseif ismember(fileName, Bill_names)
        category = 'Bill';
    elseif ismember(fileName, Nigel_names)
        category = 'Nigel';
    else
        category = 'Unknown';
    end
end
