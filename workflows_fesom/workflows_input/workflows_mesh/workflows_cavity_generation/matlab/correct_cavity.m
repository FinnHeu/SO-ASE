clear all
close all

addpath(genpath('./helpers/'))

% correct cavity mesh
meshpath='/albedo/home/fheukamp/mesh/fmesh/mesh_jigsaw_test/mesh_RTopo2.0.4_30sec_SOCAV_v5/04_correct_cavity/'
finalmeshpath='/albedo/home/fheukamp/mesh/fmesh/mesh_jigsaw_test/mesh_RTopo2.0.4_30sec_SOCAV_v5/04_correct_cavity/corrected/'

% Plotting:
plotting=0;

% Minimum number of layers:
minlayers=4;    

% Loading mesh:
fid=fopen([meshpath,'nod2d.out'],'r');
n2d=fscanf(fid,'%g',1)
nodes=fscanf(fid, '%g', [4,n2d]);
index=nodes(1,:);
xsur=nodes(2,:)';
ysur=nodes(3,:)';
nodind=nodes(4,:)';
fclose(fid);

fid=fopen([meshpath,'elem2d.out']);
el2d=fscanf(fid,'%g',1);
elem=fscanf(fid,'%g',[3 el2d]);
fclose(fid);


A=load([meshpath,'aux3d.out']);
numlevels=A(1);
Levels=A(2:numlevels+1);

depth=load([meshpath,'depth.out']);

cavity_depth=load([meshpath,'cavity_depth@node.out']);
cavity_flag=cavity_depth;
cavity_flag(cavity_flag<0)=1;

clear aux
for i=1:el2d
    aux(i)=max(abs(diff(xsur(elem(:,i)))));
end
pind=find(aux<270);
pelem=elem(:,pind);
xc=xsur(pelem);
yc=ysur(pelem);



if plotting
    figure
    tc=cavity_flag(pelem);
    pp=patch(xc,yc,tc);
    caxis([ 0 1])
    cb=colorbar;
    colormap jet
    set(pp,'EdgeColor','none');
    set(gca,'xlim',[-80 -50],'ylim',[75 83])
    set(gca,'xlim',[-90 -10],'ylim',[-85 -70])
    saveas(gcf, 'myPlot.png')
end


disp('Fixing mesh topology...')
[xsur2,ysur2,elem2,nodind2,ActiveNodes]=fix_mesh_topology(xsur,ysur,elem',logical(nodind));
n2d2=length(xsur2);
el2d2=size(elem2,1);
index2=index(ActiveNodes);
depth2=depth(ActiveNodes);
cavity_depth2=cavity_depth(ActiveNodes);
cavity_flag2=cavity_depth2;
cavity_flag2(cavity_flag2<0)=1;

xsur=xsur2;
ysur=ysur2;
elem=elem2';
nodind=nodind2;
n2d=n2d2;
el2d=el2d2;
index=index2;
depth=depth2;
cavity_depth=cavity_depth2;
cavity_flag=cavity_flag2;


clear aux
for i=1:el2d
    aux(i)=max(abs(diff(xsur(elem(:,i)))));
end
pind=find(aux<270);
pelem=elem(:,pind);
xc=xsur(pelem);
yc=ysur(pelem);

elem=elem';

disp('Correct active layers...')
FirstActiveLayer=zeros(n2d,1);
LastActiveLayer=zeros(n2d,1);
% Loop through nodes:
for ii=1:n2d
    if cavity_flag(ii)==1
    % Check if bed depth is below lowest level:
        if depth(ii)<Levels(end)
            % Use last possible layer:
            LastActiveLayer(ii)=numlevels-1; 
        else
            % Compute last active layer:
            LastActiveLayer(ii)=find(Levels<=depth(ii),1,'first')-1;
        end
        % Compute first active layer:
        FirstActiveLayer(ii)=find(Levels>=cavity_depth(ii),1,'last');
    end
end

% make deepth deeper so that we have at least 3 levels
num_layers=LastActiveLayer-FirstActiveLayer;
ind=find(num_layers<3 & num_layers>0);
depth_new=depth;
depth_new(ind)=Levels(FirstActiveLayer(ind)+3);



Depth2d=depth_new;
IceBase2d=cavity_depth;

% Compute the number of layers active in each 2D node:
% This calculation considers "layers" as depth intervals in between the
% discrete depth levels given in the Levels variable (ie, there are
% numlevels-1 layers).  A layer is active if any part of its depth range is
% above the bed or below the ice base.
% Pre-allocate:
FirstActiveLayer=zeros(n2d,1);
LastActiveLayer=zeros(n2d,1);
% Loop through nodes:
for ii=1:n2d
    % Check if bed depth is below lowest level:
    if Depth2d(ii)<Levels(end)
        % Use last possible layer:
        LastActiveLayer(ii)=numlevels-1; 
    else
        % Compute last active layer:
        LastActiveLayer(ii)=find(Levels<=Depth2d(ii),1,'first')-1;
    end
    % Compute first active layer:
    FirstActiveLayer(ii)=find(Levels>=IceBase2d(ii),1,'last');
end

%%% Ensure Overlapping Depth Levels:

Elements2d_truncated=elem;
numelements2d_truncated=length(Elements2d_truncated);
Lon2d_truncated=xsur;
Lat2d_truncated=ysur;
IsBoundary2d_truncated=nodind;
FirstActiveLayer_truncated=FirstActiveLayer;
LastActiveLayer_truncated=LastActiveLayer;
Depth2d_truncated=Depth2d;
numvertices2d_truncated=n2d;
numvertices2d=n2d;
IceBase2d_truncated=IceBase2d;

% Communicate:
disp('Deepening bed to ensure all vertex pairs have overlapping depth levels...')
t2=tic;

% Compute edge list:
Edges=edges(triangulation(Elements2d_truncated,Lon2d_truncated,Lat2d_truncated));

% Identify edges with overlapping layers:
HasOverlap=FirstActiveLayer_truncated(Edges(:,1))<=LastActiveLayer_truncated(Edges(:,2))&FirstActiveLayer_truncated(Edges(:,2))<=LastActiveLayer_truncated(Edges(:,1));

% Loop through edges:
for ii=1:size(Edges,1)
    % Check if there is overlap:
    if HasOverlap(ii)
        continue
    end
    % Identify deeper and shallower vertices:
    shallowvertexind=find(Depth2d_truncated(Edges(ii,:))==max(Depth2d_truncated(Edges(ii,:))),1,'first');
    deepvertexind=find(Depth2d_truncated(Edges(ii,:))==min(Depth2d_truncated(Edges(ii,:))),1,'first');
    % Deepen bed on shallower vertex:
    Depth2d_truncated(Edges(ii,shallowvertexind))=min(Depth2d_truncated(Edges(ii,shallowvertexind)),Levels(FirstActiveLayer_truncated(Edges(ii,deepvertexind))+1)); 
    % (min command needed in case we are deepening the same vertex
    % twice; second part of the command selected the bottom of the
    % first layer in the deeper vertex)
end
% Verify that every vertex pair now has overlapping levels:
% Pre-allocate:
FirstActiveLayer_truncated=zeros(numvertices2d,1);
LastActiveLayer_truncated=zeros(numvertices2d,1);
% Loop through nodes:
for ii=1:numvertices2d_truncated
    % Check if bed depth is below lowest level:
    if Depth2d_truncated(ii)<Levels(end)
        % Use last possible layer:
        LastActiveLayer_truncated(ii)=numlevels-1; 
    else
        % Compute last active layer:
        LastActiveLayer_truncated(ii)=find(Levels<=Depth2d_truncated(ii),1,'first')-1;
    end
    % Compute first active layer:
    FirstActiveLayer_truncated(ii)=find(Levels>=IceBase2d_truncated(ii),1,'last');
end
% Identify edges with overlapping layers:
HasOverlap=FirstActiveLayer_truncated(Edges(:,1))<=LastActiveLayer_truncated(Edges(:,2))&FirstActiveLayer_truncated(Edges(:,2))<=LastActiveLayer_truncated(Edges(:,1));
% Check overlap:
if sum(HasOverlap)<length(HasOverlap)
    error('Unable to ensure all vertex pairs have overlapping depth levels')
end


if plotting
        
    figure
    tc=depth(pelem);
    pp=patch(xc,yc,tc);
    caxis([ -800 0])
    cb=colorbar;
    colormap jet
    set(pp,'EdgeColor','none');
    set(gca,'xlim',[-80 -50],'ylim',[75 83])
    set(gca,'xlim',[-90 -10],'ylim',[-85 -70])
    saveas2('plots/mesh_so_ralph_depth_new.png',200)

    figure
    tc=depth(pelem)-Depth2d_truncated(pelem);
    pp=patch(xc,yc,tc);
    caxis([ -50 50])
    cb=colorbar;
    colormap jet
    set(pp,'EdgeColor','none');
    set(gca,'xlim',[-80 -50],'ylim',[75 83])
    set(gca,'xlim',[-90 -10],'ylim',[-85 -70])
    %saveas2('plots/mesh_so_ralph_depth_corrected_diff.png',200)

    figure
    tc=cavity_depth(pelem);
    pp=patch(xc,yc,tc);
    caxis([ -500 0])
    cb=colorbar;
    colormap jet
    set(pp,'EdgeColor','none');
    set(gca,'xlim',[-80 -50],'ylim',[75 83])
    set(gca,'xlim',[-90 -10],'ylim',[-85 -70])
    %saveas2('plots/mesh_so_ralph_cavity_depth_new.png',200)


    figure
    tc=cavity_flag(pelem);
    pp=patch(xc,yc,tc);
    caxis([ 0 1])
    cb=colorbar;
    colormap jet
    set(pp,'EdgeColor','c');
    set(gca,'xlim',[-80 -50],'ylim',[75 83])
    set(gca,'xlim',[-90 -10],'ylim',[-85 -70])
    %saveas2('plots/mesh_so_ralph_cavity_flag_new.png',200)
end

% return

% Communicate:
disp('Saving final truncated mesh...')
t3=tic;
% Write nod2d file:
disp('Writing nod2d.out...')
fid=fopen([finalmeshpath,'nod2d.out'],'w');
fprintf(fid,'%9i \n',numvertices2d);
for ii=1:numvertices2d_truncated
    fprintf(fid,'%9i %9.4f %9.4f %3i\n',ii,Lon2d_truncated(ii),Lat2d_truncated(ii),IsBoundary2d_truncated(ii));
end
fclose(fid);


% Write elem2d file:
disp('Writing elem2d.out...')
fid=fopen([finalmeshpath,'elem2d.out'],'w');
fprintf(fid,'%u \n',numelements2d_truncated);
for ii=1:numelements2d_truncated
    fprintf(fid,'%u %u %u \n',Elements2d_truncated(ii,1),Elements2d_truncated(ii,2),Elements2d_truncated(ii,3));
end
fclose(fid);

% Write aux3d file:
disp('Writing aux3d.out...')
fid=fopen([finalmeshpath,'aux3d.out'],'w');
fprintf(fid,'%8f \n',length(Levels));
for ii=1:numlevels
    fprintf(fid,'%8f \n',Levels(ii));
end
for ii=1:numvertices2d_truncated
    fprintf(fid,'%8f \n',Depth2d_truncated(ii));
end
fclose(fid);

% Write cavity depth file:
disp('Writing cavity_depth@node.out...')
fid=fopen([finalmeshpath,'cavity_depth@node.out'],'w');
for ii=1:numvertices2d_truncated
    fprintf(fid,'%8f \n',IceBase2d_truncated(ii));
end
fclose(fid);

% Save matfile:
%save([finalmeshpath,'TruncatedMesh.mat'],'*_truncated','Levels_ocean')
