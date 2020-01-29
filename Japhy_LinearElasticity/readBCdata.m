function [DBCVal,NBCVal,RBCVal] = readBCdata(idd,idn,idr)

BCfile = input('\nEnter the name of the boundary data input file(without the .inp extension): ','s');

tic;        
BCdata = fopen(strcat(BCfile,'.inp'), 'r');            
fgets(BCdata);  % Read the header line of the file
fgets(BCdata);  % Read the next dummy line.
            
% Read Dirichlet boundary data
nddf = length(find(idd == -1));  % number of dofs that specified with DBC
DBC  = zeros(nddf,3);
fgets(BCdata);  % Read the next dummy line.
for i = 1:nddf
    DBC(i,:) = str2num(fgets(BCdata));
end 
            
% Read Neumann boudnary data
nndf = length(find(idn == -1));  % number of dofs that specified with DBC
NBC = zeros(nndf,3); 
fgets(BCdata);  % Read the next dummy line.
fgets(BCdata);  % Read the next dummy line.
for i = 1:nndf
    NBC(i,:) = str2num(fgets(BCdata));
end 
            
% Read Robin boundary data
nrdf = length(find(idr == -1));  % number of dofs that specified with DBC
RBC = zeros(nrdf,4);
fgets(BCdata);  % Read the next dummy line.
fgets(BCdata);  % Read the next dummy line.
for i = 1:nrdf
    RBC(i,:) = str2num(fgets(BCdata));
end   
fclose(BCdata);
            
% Gegerate boundary value matrix
DBCVal = zeros(size(DBC,1),2);
m =1:size(DBC,1);
DBCVal(m,1) = 2*(DBC(m,1)-1) + DBC(m,2);
DBCVal(m,2) = DBC(m,3);

NBCVal = zeros(size(NBC,1),2);
m =1:size(NBC,1);
NBCVal(m,1) = 2*(NBC(m,1)-1) + NBC(m,2);
NBCVal(m,2) = NBC(m,3);

RBCVal = zeros(size(RBC,1),3);
m =1:size(RBC,1);
RBCVal(m,1) = 2*(RBC(m,1)-1) + RBC(m,2);
RBCVal(m,2) = RBC(m,3);
RBCVal(m,3) = RBC(m,4);

end   % end of function readBCdata