f=fopen('/home/loweka/testSlurmOut.txt','w+');
fprintf(f,sprintf('Current directory is %s\n',cd));
fclose('all');