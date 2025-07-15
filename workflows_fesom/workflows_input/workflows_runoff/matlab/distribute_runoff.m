% Rewrite interannual monthly runoff from Dai/Trenberth to model grid with
% the NC format.
% Add time-invariant Antarctic runoff (0.073 Sv, distributed
% along the coastline). In the data directly from Dai/Trenberth
% website, it is not included; in the data from CORE2 website, it is
% already included. We use the data from CORE2 website.
%
% Q. Wang, 06.02.2012
% adapted for JRA55, Claudia
%------------------------------------------------------------------------

clear all
close all
clc

%% ------------------------------------------------------------
% specification:
meshpath='/albedo/home/fheukamp/mesh/fmesh/mesh_jigsaw_test/mesh_RTopo2.0.4_30sec/';
meshid='DARS2.0';
datapath_river='/albedo/pool/FESOM/forcing/JRA55-do-v1.4.0_monthly_runoff/';
destpath = '/albedo/work/user/fheukamp/PostDoc2/DARS2.0/input/';

first_year=1958;
num_year=42;

extend_layer=1;   % 1: extend runoff to the specified range; 0: no extention
distribution_radius=200;   %80   %unit: km; 400 300 100 % needs to be larger than the largest element

r_earth=6367.5e3;    %m, FEOM fixed parameter

%% ---------------------------------------------------------
% Loading mesh:
loadmesh=1;
if loadmesh
    disp('load mesh')
    
    fid=fopen([meshpath,'nod2d.out']);
    n2d=fscanf(fid,'%g',1);
    nodes=fscanf(fid, '%g', [4,n2d]);
    fclose(fid);
    
    xcoord=nodes(2,:);
    ycoord=nodes(3,:);
    nodind=nodes(4,:);
    
    fid=fopen([meshpath,'elem2d.out']);
    el2d=fscanf(fid,'%g',1);
    elem=fscanf(fid,'%g',[3 el2d]);
    fclose(fid);
end

coast = find(nodind==1);

deg2rad=pi/180;
xc = xcoord(coast)*deg2rad;   % only coastline nodes
yc = ycoord(coast)*deg2rad;

%-----------------------------------------------------------
% fesom cluster arear

disp('calculate cluster area')

domain_length=2*pi;  %360*deg2rad;

patch_area=zeros(1,n2d);
%h = waitbar(0, 'Please wait... Looping through elements.');
for i=1:el2d;

 %   waitbar(i/el2d, h)
    elnode=elem(:,i);
    local_cart(1,:)=xcoord(elnode)*deg2rad;
    local_cart(2,:)=ycoord(elnode)*deg2rad;
    
    %  T R A N S F O R M A T I O N - M A T R I X   "jacobian"
    jacobian2D = local_cart(:,[2 3])-local_cart(:,[1 1]);
    
    % check cyclic boundary
    for j=1:2
        if (jacobian2D(1,j)> domain_length/3.0)
            jacobian2D(1,j)=jacobian2D(1,j)-domain_length; end;
        if (jacobian2D(1,j)<-domain_length/3.0)
            jacobian2D(1,j)=jacobian2D(1,j)+domain_length; end;
    end
    
    jacobian2D=jacobian2D*r_earth;
    %midlat=mean(local_cart(2,:));
    jacobian2D(1,:)=jacobian2D(1,:)*mean(cos(local_cart(2,:)));
    determinant=det(jacobian2D);
    voltriangle = abs(determinant)/2.0;
    
    patch_area(elnode)=patch_area(elnode)+voltriangle;

    % Print progress every 5% (or adjust as needed)
    if mod(i, round(el2d/100)) == 0 || i == el2d
        progress = (i / el2d) * 100;
        fprintf('Progress: %.0f%%\n', progress);
    end
end
%delete(h)
patch_area=patch_area/3;

% Build nod_in_elem2D and nghbr_nod2d



deg2rad=pi/180;
rearth=6367.5;  %km
if(extend_layer)
    disp('build nghbr_nod2d')
    
    h = waitbar(0, 'Please wait... Looing through nods.');
    for node=1:n2d
        waitbar(node/n2d,h)
        if mod(node, 1e3) == 0
            progress = (node * 100) / n2d;
            fprintf('Progress: %.0f%%\n', progress);
            %disp(num2str(node*100/n2d),'%')
        end
        % to spped up, only build for coastal nodes
        if(nodind(node)==0)
            continue,
        end
        
        ln = deg2rad*xcoord;
        lt = deg2rad*ycoord;
        xn=ln(node);
        yn=lt(node);
        
        alpha = ...
            acos( cos(yn)*cos(lt).*cos(xn-ln) ...
            + sin(yn)*sin(lt) );
        dx = rearth*abs(alpha);
        
                
        xn=xcoord(node);
        yn=ycoord(node);
        
        %if(xn>-51.5 & xn<-48 & yn>-1 & yn<2) 
        %    a=find(dx<=distribution_radius*2);  %km
        %else
            a=find(dx<=distribution_radius);  %km
        %end
            
        nghbr_nod2d(node).addresses=a;
        dist_nod2d(node).addresses=dx(a);

        
    end
    
end

%--------------------------------------------------------

% cell area and coordinates on the source data grid

ni=1440;
nj=720;


longitude = ncread([datapath_river,'runoff_cell_area.15Dec2016.nc'],'longitude'); %degree
latitude = ncread([datapath_river,'runoff_cell_area.15Dec2016.nc'],'latitude'); %degree
area = ncread([datapath_river,'runoff_cell_area.15Dec2016.nc'],'areacello');% m^2

[xreg,yreg]=meshgrid(longitude,latitude);
xreg=xreg';
yreg=yreg';

xreg=xreg*deg2rad;  % cell center coordinates to radian
yreg=yreg*deg2rad;


%% ---------------------------------------------------------

river_runoff_fesom = zeros(1,n2d);
river_runoff_coast = zeros(1,length(coast));
aux_runoff=river_runoff_fesom;

for year=first_year:first_year+num_year-1    
    
    disp(num2str(year))
    
    %--------------------------------
    riverfile=[datapath_river,'friver_monthly.',num2str(year),'.nc'];
    
    time=ncread(riverfile,'time');
    num_days=length(time);
    %fmesh
    % new output file
    
    % one file per year
    %fid=netcdf.create([meshpath,'JRA55_runoff/', num2str(year), '.nc'],'NC_CLOBBER');
    fid=netcdf.create([destpath,'forcing_data_on_grid/', num2str(year), '.nc'],'NC_CLOBBER');
    
    % define variables
    daydimID = netcdf.defDim(fid,'Dimday',num_days);
    nod2ddimID = netcdf.defDim(fid,'Dimnod2d',n2d);
    rr_ID = netcdf.defVar(fid,'runoff','double',[nod2ddimID, daydimID]);
    % time at the end
    netcdf.endDef(fid);
    

    
    %varid_runoff = netcdf.inqVarID(ncid,'runoff');
    for day=1:num_days
        
        day
        % reading river runoff data
        
        start=[1 1 day];
        count=[ni nj 1];
        stride=[1 1 1];
        runoff_raw = ncread(riverfile,'friver',start,count,stride); % kg/m2/s
        
        % convert to kg/s
        runoff_raw=runoff_raw.*area;
        
        %--------------------------------------
        
        river_runoff_coast(:)=0.0;
        river_runoff_fesom(:)=0.0;
        runoff_cell=find(runoff_raw>0);
        
        for k=1:length(runoff_cell)
            
            % coordinates of river runoff boxes (centered)
            yk=yreg(runoff_cell(k));
            xk=xreg(runoff_cell(k));
            
            % find next coastal grid node
            dist = ...
                acos( cos(yk)*cos(yc).*cos(xk-xc) ...
                + sin(yk)*sin(yc) );
            %dist = r_earth*abs(alpha);
            dist=abs(dist);
            
            [minval, node]=min(dist);
            
            river_runoff_coast(node) = river_runoff_coast(node) +  ...
                runoff_raw(runoff_cell(k));
        end
        
        % kg/m2/s
        river_runoff_coast = river_runoff_coast./patch_area(coast);
        
        % put to fesom grid
        river_runoff_fesom(coast)=river_runoff_coast;
        
        % check
        disp([num2str(sum(river_runoff_fesom.*patch_area)/1e9)]);
        disp([num2str(sum(sum(runoff_raw))/1e9)]);
        
        %--------------------------------------
        
        if(extend_layer)
            %disp('extend river mouth ...')
            aux_runoff=river_runoff_fesom;
            river_runoff_fesom(:)=0;
            
            for node=1:n2d
                if(aux_runoff(node)==0) continue, end
                               
                
                dmax=max(dist_nod2d(node).addresses);
                weight=patch_area(nghbr_nod2d(node).addresses).* (1-dist_nod2d(node).addresses/dmax);
                %weight=patch_area(nghbr_nod2d(node).addresses).* abs((1-dist_nod2d(node).addresses/distribution_radius));
                
                if weight==0
                    weight
                end
                weight=weight/sum(weight);
                
                river_runoff_fesom(nghbr_nod2d(node).addresses) = ...
                    river_runoff_fesom(nghbr_nod2d(node).addresses) + ...
                    aux_runoff(node)*patch_area(node).*weight;
            end
            
            river_runoff_fesom=river_runoff_fesom./patch_area;
            
            disp([num2str(sum(river_runoff_fesom.*patch_area)/1e9)])
            
        end
        
        %--------------------------------------
        
        % save to output files
        
        % put data
        istart=[0,   day-1];
        icount=[n2d, 1];
        netcdf.putVar(fid, rr_ID, istart, icount, river_runoff_fesom);
        
    end     % month
    
    % close file
    netcdf.close(fid);
    
end     %year


%-------------------------------------------------------------------------

return
