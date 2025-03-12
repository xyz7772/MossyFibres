function category = getMice(fileName)
    Initialize
    if ismember(fileName, A1_data)
        category = 'Animal1';
    elseif ismember(fileName, A2_data)
        category = 'Animal2';
    elseif ismember(fileName, A3_data)
        category = 'Animal3';
    elseif ismember(fileName, A4_data)
        category = 'Animal4';
    else
        category = 'Unknown';
    end
end
