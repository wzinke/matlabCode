function klParSavev1(fName,varVal,varName),
    eval(sprintf('%s=varVal;',varName));
    save(fName,varName,'-v7.3');
end