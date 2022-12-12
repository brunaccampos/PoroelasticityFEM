function WriteMesh2VTK(filename, description, Mesh, scalardata, vectordata)
% Write FEM mesh data to a file in VTK format
% ------------------------------------------------------------------------
% Inputs
% ------------------------------------------------------------------------
% filename = string, should end with .vtk
% description = string describing the model
% Mesh = displacement or pressure field
% scalardata = data structure
% ------------------------------------------------------------------------
% Creator: Robert Gracie, University of Waterloo 2012, 2018.
% Use of this file is governed under the GNU Public Licence:
% https://www.gnu.org/licenses/gpl-3.0.en.html
% ------------------------------------------------------------------------
% Last modified: June 2022
% Modified by: Bruna Campos
% Adapt for poroelasticity model

%% Mesh properties
nn = Mesh.nn; % number of nodes
nsd = Mesh.nsd; % number of spatial dimensions
ne = Mesh.ne; % number of elements
nne = Mesh.nne; % number of nodes per element
nd = length(scalardata); % number of scalar datum defined at each node
nv = length(vectordata); % number of vector datum defined at each node
nodes = Mesh.coords; % nodal coordinates
conn = Mesh.conn; % mesh connectivity

% fix length of nodal coordinates vector
if nsd == 1 && size(nodes,2) == 1
    nodes = [nodes, zeros(nn,1), zeros(nn,1)];
elseif nsd == 2 && size(nodes,2) == 2
    nodes = [nodes, zeros(nn,1)];
end

min_sig_digs = 5;
max_sig_digs = 20;

%% File header
% open the file in write mode
fid = fopen(filename, 'w');
% header
fprintf(fid,'%s\n','# vtk DataFile Version 2.0');
fprintf(fid,'%s\n',['HFX mesh description: ',description]);
fprintf(fid,'%s\n','ASCII');
fprintf(fid,'%s\n','DATASET UNSTRUCTURED_GRID');

%% Element connectivity
if nne == 2 % L2
    outputformat = '%d %d %d  \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 3;
elseif nne == 3 && nsd == 1 % L3
    outputformat = '%d %d %d %d \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 21;
elseif nne == 4 % Q4 FOR PRESSURE FIELD
    outputformat = '%d %d %d %d %d \n';
    %   conn = [4*ones(ne,1),conn(:,1:4)-ones(ne,4)];
    conn = [4*ones(ne,1),conn-ones(ne,4)];
    cell_type = 9;
elseif nne == 9 % Q9
    outputformat = '%d %d %d %d %d %d %d %d %d \n';
    middlenodes = conn(:,9);
    conn(:,9) = [];
    % remove middle nodes from node list
    nn = nn - length(middlenodes);
    nne = 8;
    nodes(middlenodes,:) = [];
    % remove middle nodes from scalar data
    ndata = length(scalardata);
    for field = 1:ndata
        fielddata = scalardata(field).data;
        fielddata(middlenodes,:) = [];
        scalardata(field).data = fielddata;
    end
    % remove middle nodes from vector data
    if ~isempty(vectordata)
        ndata = length(vectordata);
        for field = 1:ndata
            fielddata = vectordata(field).data;
            fielddata(middlenodes,:) = [];
            vectordata(field).data = fielddata;
        end
    end   
    % remove middle nodes from connectivity
    middlenodes = sort(middlenodes,'ascend');
    for i = 1:length(middlenodes)
        conn(conn>middlenodes(i)) = conn(conn>middlenodes(i))-1;
        middlenodes(middlenodes>middlenodes(i)) = middlenodes(middlenodes>middlenodes(i))-1;
    end   
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 23;
    elseif nne == 3 && nsd == 2 % T3
    outputformat = '%d %d %d %d \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 5;
    elseif nne == 6 % T6
    outputformat = '%d %d %d %d %d %d \n';
    conn = [nne*ones(ne,1),conn-ones(ne,nne)];
    cell_type = 22;
end

%% Print data
% nodal coordinates
text = ['POINTS ',num2str(nn),' float'];
fprintf(fid,'%s\n',text);
fprintf(fid,'%f %f %f\n',nodes');
% cells: # elements, # nodes
text = ['CELLS ',num2str(ne), ' ',num2str((nne+1)*ne)];
fprintf(fid,'%s\n',text);
% connectivity
fprintf(fid,outputformat,conn');
% cell types
text = ['CELL_TYPES ',num2str(ne)];
fprintf(fid,'%s\n',text);
% number of nodes per cell for each element
fprintf(fid,'%d\n',cell_type*ones(1,ne));

% checking whether point data is defined
if nd > 0
    % print POINT_DATA and # of nodes for each data set
    text = ['POINT_DATA ',num2str(nn)];
    fprintf(fid,'%s\n',text);
end

% loop over scalar data sets
for i = 1:nd
    if isfield(scalardata,'type')
        type = scalardata(i).type;
    else
        type = 'float';
    end
    
    if strcmp(type,'float')
        sig_digs = real(floor(log10(scalardata(i).data)));
        lowest_exp = min(sig_digs(~isinf(sig_digs)));
        numsigdigs = min(max_sig_digs,max([abs(lowest_exp),min_sig_digs]));
    elseif strcmp(type,'int')
        numsigdigs = '0';
    end
    % number of components
    numcomp = size(scalardata(i).data,2);
    
    if numcomp == 3
        text = ['VECTORS ', scalardata(i).name,' ' type ' '];
        fprintf(fid,'%s\n',text);
    else
        text = ['SCALARS ',scalardata(i).name,' ' type ' ', num2str(numcomp)];
        fprintf(fid,'%s\n',text);
        text = 'LOOKUP_TABLE default';
        fprintf(fid,'%s\n',text);
    end
    fprintf(fid,[repmat(['%.' num2str(numsigdigs) 'f '],1,numcomp),'\n'],scalardata(i).data');
end

% loop over vector data sets
if nv > 0
    for i = 1:nv
        if isfield(vectordata,'type')
            type = vectordata(i).type;
        else
            type = 'float';
        end
        
        if strcmp(type,'float')
            sig_digs = real(floor(log10(vectordata(i).data)));
            lowest_exp = min(sig_digs(~isinf(sig_digs)));
            numsigdigs = min(max_sig_digs,max([abs(lowest_exp),min_sig_digs]));
        elseif strcmp(type,'int')
            numsigdigs = '0';
        end
        
        text = ['VECTORS ', vectordata(i).name,' ' type ' '];
        fprintf(fid,'%s\n',text);
        for j = 1:nn
            fprintf(fid,[repmat(['%.' num2str(numsigdigs) 'f '],1,numcomp),'\n'],vectordata(i).data(j,1),vectordata(i).data(j,2), 0);
        end
    end
end

% close file
fclose(fid);

end