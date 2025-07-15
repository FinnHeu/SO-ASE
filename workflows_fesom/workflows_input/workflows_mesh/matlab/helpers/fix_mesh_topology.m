function [X2d,Y2d,Elements2d,IsBoundary2d,ActiveNodes]=fix_mesh_topology(X2d,Y2d,Elements2d,IsBoundary2d,varargin)

% FixMeshTopology
% Mike Wolovick, 3/25/2024

% This function fixes the topology of an input mesh by removing vertices
% and elements where necessary.

% "Fixing the mesh topology" means removing bottleneck nodes, edges with
% two boundary nodes (ie, bottleneck elements), isolated lakes, and then 
% iteratively removing elements with only one neighbor.
% Note:  by removing bottleneck nodes before checking for lakes, it is
% not necessary to do a separate check for isolated inlets.

% The final output mesh produced by this function is a subset of the
% original input mesh that is well-connected, without single-node or
% single-element bottlenecks.

% Note that the "ActiveNodes" output is the same size as the original
% vertex list.  It can thus be used to extract information that was stored
% on the original vertices other than coordinates.  For example:
% Variable_old=[some information on the old vertex list];
% [X2d,Y2d,Elements2d,IsBoundary2d,ActiveNodes]=FixMeshTopology(X2d,Y2d,Elements2d,IsBoundary2d)
% Variable_new=Variable_old(ActiveNodes);

% This function calls TruncateMesh to handle the indexing of the new
% elements list when removing vertices or elements from a mesh.

% Optional input arguments: (name-value pairs)
% 'cutbottlenecknodes'      logical
% 'cutbottleneckelements'   logical
% 'cutlakes'                logical
% 'cutoneneighborelements'  logical

%% Arbitrary Parameters:

% This function contains a few arbitrary parameters not set by the input.
% Those are:
relmaxnumconnections=20; %10           % maximum number of connections to try when searching for contiguous regions, in multiples of sqrt(numvertices2d)
maxnumseedpoints=20;               % maximum number of random seed points to consider when searching for the main mesh region (as opposed to the lakes)
numsearchseedpoints=10;            % number of random seed points to consider before taking the biggest contiguous mesh region we can find
maxoneneighboriterations=100;      % maximum number of iterations when removing elements with only one neighbor

%% Input Checking:

% Check input vertices lists:
% check x-coordinate size:
if size(X2d,1)>1 && size(X2d,2)>1
    error('Input x-coordinate list must be a vector')
elseif size(X2d,2)>1
    X2d=X2d';
end
% Compute number of vertices:
numvertices2d=length(X2d);
% Check size of other vertices lists:
if size(Y2d,1)>1 && size(Y2d,2)>1
    error('Input x-coordinate list must be a vector')
elseif size(Y2d,2)>1
    Y2d=Y2d';
end
if length(Y2d)~=numvertices2d
    error('Input x-corrdinate and y-coordinate lists must be the same size')
end
if size(IsBoundary2d,1)>1 && size(IsBoundary2d,2)>1
    error('Input IsBoundary2d list must be a vector')
elseif size(IsBoundary2d,2)>1
    IsBoundary2d=IsBoundary2d';
end
if length(IsBoundary2d)~=numvertices2d
    error('Input IsBoundary2d list must be the same length as the coordiante lists')
end
% Check contents of vertices lists:
if sum(isnan(X2d))>0 || sum(isinf(X2d))>0 
    error('Input x-coordinate list cannot contain NaN or Inf')
elseif sum(isnan(Y2d))>0 || sum(isinf(Y2d))>0 
    error('Input y-coordinate list cannot contain NaN or Inf')
end
if ~islogical(IsBoundary2d)
    error('Input IsBoundary2d must be type "logical".')
end

% Check element indexing:
% Check size of elements list and compute number of elements:
if size(Elements2d,2)~=3
    if size(Elements2d,1)==3
        Elements2d=Elements2d';
        numelements2d=size(Elements2d,1);
    else
        error('Input element indexing must be either a 3-column or a 3-row matrix')
    end
else
    numelements2d=size(Elements2d,1);
end
% Check contents of elements list:
if sum(isnan(Elements2d(:)))>0  || sum(isinf(Elements2d(:)))>0 
    error('Input element indexing cannot contain NaN or Inf')
elseif max(Elements2d(:))>numvertices2d
    error('Input element indexing cannot contain values larger than the number of input vertices')
elseif min(Elements2d(:))<1
    error('Input element indexing cannot contain values less than 1')
elseif ~isequal(Elements2d,round(Elements2d))
    error('Input element indexing must contain only integers')
end

%% Parse Optional Inputs:

% Assign default values:
cutbottlenecknodes=1;
cutbottleneckelements=1;
cutlakes=1;
cutoneneighborelements=1;

% Check for additional input arguments:
if ~isempty(varargin)
    % Check length of additional arguments:
    if length(varargin)~=2 && length(varargin)~=4 && length(varargin)~=6 && length(varargin)~=8 
        error('Incorrect number of input arguments')
    end
    % Parse name/value pairs:
    for ii=1:length(varargin)/2
        if strcmp(varargin{(ii-1)*2+1},'cutbottlenecknodes')
            cutbottlenecknodes=varargin{ii*2};
        elseif strcmp(varargin{(ii-1)*2+1},'cutbottleneckelements')
            cutbottleneckelements=varargin{ii*2};
        elseif strcmp(varargin{(ii-1)*2+1},'cutlakes')
            cutlakes=varargin{ii*2};
        elseif strcmp(varargin{(ii-1)*2+1},'cutoneneighborelements')
            cutoneneighborelements=varargin{ii*2};
        else
            if ischar(varargin{(ii-1)*2+1})
                error(['Unrecognized input parameter, "',varargin{(ii-1)*2+1},'"',', recognized values are "cutbottlenecknodes", "cutbottleneckelements", "cutlakes", and "cutoneneighborelements".'])
            else
                error('Unable to parse additional input parameters')
            end
        end
    end
end

%% Preparation:

% Compute maximum number of connections:
maxnumconnections=ceil(sqrt(numvertices2d)*relmaxnumconnections);

% Get all necessary mesh connectivity information:  
% Make a triangulation object:
ThisTriangulation=triangulation(Elements2d,X2d,Y2d);
% Get list of elements attached to each node:
AttachedElements=vertexAttachments(ThisTriangulation);
% Generate list of element neightbors:
ElementNeighbors=neighbors(ThisTriangulation);

% Note: the above procedure lists each node as connected to itself, but
% that is fine for my purposes.

%% Identify and Process Bottlenecks:

% Identify bottleneck nodes:
BottleneckNodes=false(numvertices2d,1);
% Check whether to cut bottleneck nodes:
if cutbottlenecknodes
    % Loop through starting vertices:
    for ii=1:numvertices2d
        % Check if this is a boundary node with more than 1 attached element:
        if ~IsBoundary2d(ii) || length(AttachedElements{ii})==1
            continue
        end
        % Loop through attached elements to look for odd-one-out bottleneck:
        for jj=1:length(AttachedElements{ii})
            % Check if this element is bottlenecked at the original node:
            % (ie, look for an element that is not a neighbor to any of the
            % other elements attached to this node)
            % NOTE: this bottleneck-checking cannot detect a situation
            % where two groups of two elements are bottlenecked at the
            % node.  This code can only detect the situation where one
            % side of the bottleneck is served by a single element.
            % Check the number of unique element neighbors:
            if length(unique([ElementNeighbors(AttachedElements{ii}(jj),:),AttachedElements{ii}]))==3+length(AttachedElements{ii})  % if any elements are neighbors with each other, this list will be shorter than maximum length
                % Flag this vertex and break from inner loop:
                BottleneckNodes(ii)=1;
                break
            end
        end
        % Check for 2vs2 bottleneck:
        if length(AttachedElements{ii})==4 && ~BottleneckNodes(ii)
            % Check 3 possible 2vs2 combinations:
            possiblecombos=[1,2,3,4;1,3,2,4;1,4,2,3]; % each row is [group1,group2]
            for jj=1:size(possiblecombos,1)
                % Check if the first group has any neighbors in the second
                % group:
                list1=ElementNeighbors(AttachedElements{ii}(possiblecombos(jj,1:2)),:);
                list1=list1(:);
                list2=AttachedElements{ii}(possiblecombos(jj,3:4));
                if length(unique([list1;list2(:)]))==length(list1)+length(list2)
                    % Flag this vertex and break from inner loop:
                    BottleneckNodes(ii)=1;
                    break
                end
            end
        end
    end
end

% Check whether to cut bottleneck elements:
if cutbottleneckelements
    % Identify and deactivate bottlenecked elements:
    ActiveElements=~(sum(BottleneckNodes(Elements2d),2)>0|sum(IsBoundary2d(Elements2d),2)==3);
    % Note that an element with 3 boundary vertices must always be a bottleneck
    % between its connected elements.
else
    ActiveElements=sum(BottleneckNodes(Elements2d),2)==0;
end

%% Make a new mesh without the bottlenecks:

% Truncate the mesh
[Elements2d,ActiveNodes1,ThisIsBoundary]=TruncateMesh(Elements2d(ActiveElements,:),true(numvertices2d,1));
% Truncate the coordinate lists:
X2d=X2d(ActiveNodes1);
Y2d=Y2d(ActiveNodes1);
% Compute boundary nodes:
IsBoundary2d=ThisIsBoundary|IsBoundary2d(ActiveNodes1);
% Get new mesh size:
numvertices2d=length(X2d);
numelements2d=size(Elements2d,1);

% Get all necessary mesh connectivity information:
% Make a triangulation object:
ThisTriangulation=triangulation(Elements2d,X2d,Y2d);
% Generate list of element neightbors:
ElementNeighbors=neighbors(ThisTriangulation);
% Generate list of edges:
Edges=edges(ThisTriangulation);

%% Remove disconnected lakes:

% Check whether to cut lakes:
if cutlakes
    % Pre-allocate number of connections record:
    NumConnections=zeros(maxnumseedpoints,1);
    % Iteratively test seed points:
    done_outer=0;
    counter_outer=1;
    while done_outer==0
        % Reset active nodes variable:
        ActiveNodes2=false(numvertices2d,1);
        % Randomly select a seed point:
        ActiveNodes2(randi(numvertices2d))=1;
        % Iteratively trace connections from the seed point:
        done_inner=0;
        counter_inner=1;
        while done_inner==0
            % Remember old list of connections:
            ActiveNodes_last=ActiveNodes2;
            % Identify (new) connected edges:
            ConnectedEdges=sum(ActiveNodes2(Edges),2)==1;
            % Add new connected nodes:
            ActiveNodes2(Edges(ConnectedEdges,:))=1;
            % Check if we've made new connections and break from loop:
            if isequal(ActiveNodes2,ActiveNodes_last)
                done_inner=1;
            elseif counter_inner>maxnumconnections
                error('Unable to identify disconnected lakes (unable to trace all connections from seed point)')
            else
                counter_inner=counter_inner+1;
            end
        end
        % Record number of connections:
        NumConnections(counter_outer)=sum(ActiveNodes2);
        % Check whether to break from loop:
        if NumConnections(counter_outer)>numvertices2d/2 % if we connected to more than half, then we know that this is the main region
            done_outer=1;
        elseif NumConnections(counter_outer)==max(NumConnections(1:counter_outer)) && counter_outer>numsearchseedpoints % if we've gone more than numsearchseedpoints, then just take the biggest region we have found
            done_outer=1;
        elseif counter_outer>maxnumseedpoints
            error('Unable to identify disconnected lakes (unable to find the right seed point, mesh may have many disconnected regions)')
        else
            counter_outer=counter_outer+1;
        end
    end
    % Deactivate disconnected nodes and elements:
    ActiveElements=sum(ActiveNodes2(Elements2d),2)==3;
else
    ActiveNodes2=true(numvertices2d,1);
    ActiveElements=true(numelements2d,1);
end

%% Cut Out Elements with One or Fewer Neighbors:

% Check whether to cut one-neighbor elements:
if cutoneneighborelements
    % Pre-allocate:
    ActiveNeighbors=false(numelements2d,3);
    % Iterate:
    done=0;
    counter=1;
    while done==0
        % Identify active neighbors:
        ActiveNeighbors(~isnan(ElementNeighbors))=ActiveElements(ElementNeighbors(~isnan(ElementNeighbors)));
        % Identify elements to deactivate:
        TooFewNeighbors=sum(ActiveNeighbors,2)<=1&ActiveElements;
        % Deactivate them:
        ActiveElements(TooFewNeighbors)=0;
        % Break from loop:
        if sum(TooFewNeighbors)==0
            done=1;
        elseif counter>maxoneneighboriterations
            error('Unable to remove all elements with 1 or fewer neighbors')
        else
            counter=counter+1;
        end
    end
end

%% Make the final output mesh:

% Truncate the mesh
[Elements2d,ActiveNodes2,ThisIsBoundary]=TruncateMesh(Elements2d(ActiveElements,:),ActiveNodes2);
% Truncate the coordinate lists:
X2d=X2d(ActiveNodes2);
Y2d=Y2d(ActiveNodes2);
% Compute boundary nodes:
IsBoundary2d=ThisIsBoundary|IsBoundary2d(ActiveNodes2);
% Compute final active nodes list:
ActiveNodes=ActiveNodes1;
ActiveNodes(ActiveNodes)=ActiveNodes2;
