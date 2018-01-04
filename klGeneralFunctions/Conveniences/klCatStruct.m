function catStruct = klCatStruct(old,new,dim)

fNames = fieldnames(old);

for iField = 1:length(fNames)
    if isstruct(old.(fNames{iField}))
        catStruct.(fNames{iField}) = klCatStruct(old.(fNames{iField}),new.(fNames{iField}),dim);
    elseif ischar(old.(fNames{iField}))
        catStruct.(fNames{iField}) = old.(fNames{iField});
    else
        catStruct.(fNames{iField}) = cat(dim,old.(fNames{iField}),new.(fNames{iField}));
    end
end