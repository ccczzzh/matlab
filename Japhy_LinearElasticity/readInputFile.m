function [ndm,ndf,eltype,nel,np,nen,elpar,x0,...
                          elc,ndbc,nnbc,nrbc,vdbc,vnbc,vrbc] = readInputFile()

prName = input('\nEnter the name of the input file(without the .inp extension): ','s');
inputFile = fopen(strcat(prName,'.inp'), 'r');
            
fgets(inputFile);  % Read the header line of the file
fgets(inputFile);  % Read the next dummy line.
            
ndm      = fscanf(inputFile, 'ndm      :%d');    fgets(inputFile);
ndf     = fscanf(inputFile,  'ndf      :%d');    fgets(inputFile);
eltype   = fscanf(inputFile, 'eltype   :%d');    fgets(inputFile);
nel    = fscanf(inputFile,   'nel      :%d');    fgets(inputFile);
np     = fscanf(inputFile,   'np       :%d');    fgets(inputFile);
nen      = fscanf(inputFile, 'nen      :%d');
fgets(inputFile);  fgets(inputFile);
elpar = str2num(fgets(inputFile));
               
% Read node coordinates
x0 = zeros(np,ndm);
fgets(inputFile);  % Read the next dummy line.
fgets(inputFile);  % Read the next dummy line.
for m = 1:np
    dummy = str2num(fgets(inputFile)); 
    x0(m,1) = dummy(2);
    x0(m,2) = dummy(3);
end 
            
% Read connectivity of the elements, i.e. ix
elc = zeros(nen,nel);
fgets(inputFile);  % Read the next dummy line.
fgets(inputFile);  % Read the next dummy line.
for in = 1:nen
    dummy = str2num(fgets(inputFile)); %#ok<*ST2NM>
    elc(in,:) = dummy(2:nel+1);  
end
            
% Read the number of nodes with different boundary conditions
fgets(inputFile);  % Read the next dummy line.
ndbc = fscanf(inputFile, 'ndbcnodes :%d');  fgets(inputFile);
nnbc = fscanf(inputFile, 'nnbcnodes :%d');  fgets(inputFile);
nrbc = fscanf(inputFile, 'nrbcnodes :%d');  fgets(inputFile);
            
% Read DBC nodes and DOFs
vdbc = zeros(ndbc,ndf+1);
fgets(inputFile);  % Read the next dummy line.
fgets(inputFile);  % Read the next dummy line.
for m = 1:ndbc
    vdbc(m,:) = str2num(fgets(inputFile));
end 
            
% Read NBC nodes and DOFs
vnbc = zeros(nnbc,ndf+1);
fgets(inputFile);  % Read the next dummy line.
fgets(inputFile);  % Read the next dummy line.
for m = 1:nnbc
    vnbc(m,:) = str2num(fgets(inputFile));
end

% Read RBC nodes and DOFs
vrbc = zeros(nrbc,ndf+1);
fgets(inputFile);  % Read the next dummy line.
fgets(inputFile);  % Read the next dummy line.
for m = 1:nrbc
    vrbc(m,:) = str2num(fgets(inputFile));
end

fclose(inputFile);   % Close the input data file.

end   % end of function readInputFile