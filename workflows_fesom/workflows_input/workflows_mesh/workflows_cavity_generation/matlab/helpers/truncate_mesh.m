function [Elements2d_new,ActiveNodes_old,IsBoundary_new]=truncate_mesh(Elements2d_old,ActiveNodes_old)

% TruncateMesh
% Mike Wolovick, 3/25/2024

% This function performs the indexing operations for removing nodes from a
% mesh and updating the elements array.  Note that the final list of the 
% new mesh may not include all of the requested vertices, because this 
% function removes disconnected nodes (those that are not part of an active
% element).  Thus, an updated ActiveNodes list is part of the output.

% suffixes "_new" and "_old" refer to the new and old meshes.

% Inputs: 
% Elements2d_old      connectivity information of old mesh
% ActiveNodes_old     logical vector indicating which nodes to include in the
%                      new mesh

% Outputs:
% Elements2d_new      connectivity of new mesh
% ActiveNodes_old     revised list of active nodes (note that size matches
%                     the OLD node list)
% IsBoundary_new      list of new boundary nodes in the new mesh (note that
%                     this will not flag nodes that were alo boundaries in 
%                     the old mesh, only new boundaries) 

% In order to get a complete list of boundary nodes in the new mesh, you
% must evaluate:
% IsBoundary_new=IsBoundary_new|IsBoundary_old(ActiveNodes_old);

%% Input Checking:

% Check active nodes list:
if size(ActiveNodes_old,1)>1 && size(ActiveNodes_old,2)>1
    error('Input active nodes list must be a vector')
elseif size(ActiveNodes_old,2)>1
    ActiveNodes_old=ActiveNodes_old';
end
if ~islogical(ActiveNodes_old)
    error('Input active nodes list must be type "logical".')
end

% Compute size of old mesh:
numvertices2d_old=length(ActiveNodes_old);

% Check element indexing:
if size(Elements2d_old,2)~=3
    if size(Elements2d_old,1)==3
        Elements2d_old=Elements2d_old';
    else
        error('Input element indexing must be either a 3-column or a 3-row matrix')
    end
end
if sum(isnan(Elements2d_old(:)))>0 || sum(isinf(Elements2d_old(:)))>0 
    error('Input element indexing cannot contain NaN or Inf')
elseif max(Elements2d_old(:))>numvertices2d_old
    error('Input element indexing cannot contain values larger than the number of input vertices')
elseif min(Elements2d_old(:))<1
    error('Input element indexing cannot contain values less than 1')
elseif ~isequal(Elements2d_old,round(Elements2d_old))
    error('Input element indexing must contain only integers')
end

%% Work:

% Identify active elements:
ActiveElements=sum(ActiveNodes_old(Elements2d_old),2)==3;  % only include elements with 3 active nodes

% Identify disconnected nodes:
ConnectedNodes=false(numvertices2d_old,1);
ConnectedNodes(Elements2d_old(ActiveElements,:))=1;

% Modify active node list to remove disconnected nodes:
ActiveNodes_old=ActiveNodes_old&ConnectedNodes;

% Identify elements just outside the new mesh:
JustOutsideElements=sum(ActiveNodes_old(Elements2d_old),2)>0&sum(ActiveNodes_old(Elements2d_old),2)<3;

% Compute size of new mesh:
numvertices2d_new=sum(ActiveNodes_old);

% Identify boundary nodes:
IsBoundary_new=false(numvertices2d_old,1);  % start the same size as the old mesh
IsBoundary_new(Elements2d_old(JustOutsideElements,:))=1;
IsBoundary_new=IsBoundary_new(ActiveNodes_old);  % truncate down to new mesh size

% Generate reverse indexing key for new elements list:
% Record original indices of truncated list:
OldInd_old=linspace(1,numvertices2d_old,numvertices2d_old)';
OldInd_new=OldInd_old(ActiveNodes_old);
% Pre-allocate:
IndexingKey_old=zeros(numvertices2d_old,1);  % size of array matches old list, value indicates index in new list
% Loop through new nodes:
for ii=1:numvertices2d_new
    % Assign reverse index value:
    IndexingKey_old(OldInd_new(ii))=ii;
end

% Generate new elements list:
Elements2d_new=IndexingKey_old(Elements2d_old(ActiveElements,:));
