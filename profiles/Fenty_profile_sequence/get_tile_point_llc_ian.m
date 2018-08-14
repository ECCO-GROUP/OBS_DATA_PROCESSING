function [MITprof] = ...
    get_tile_point_llc_ian(XC, YC, ni, nj, nx, MITprof)
%
% This script finds the important tile coordinates for the mitgcm profile
% package for a profile point on an llc grid.
%
% Filename; get_tile_point_llc_ian
%  ** former filename : updates step04a_get_tile_point [from An]
% 
% Date Created: 2014-07-30
% Last Modified: 
%
% notes:
%    previous version is not called as a function but as a script
%    here, the critical dependence on prof_point is made explicit.
%
%
% function [MITprof] = ...
%   get_tile_point_llc_ian(lon_grid, lat_grid, ni, nj, nx, MITprof)
%
%     XC , YC                : the XC and YC format (not patchface)
%     ni,nj                  : the tile size for the model
%     nx                     : the llcNX size
%     MITprof                : the MIT profile structure containing the
%       important filed 'prof_point' --> the index corresponding to  each
%       profile's nearest neighbor point in the llc grid [in patchface format]
%       this prof_point is given by get_profpoint_llc_ian
%
% notes:
% ni=30;nj=30;	% tile size for v4
%
%  FROM SIZE.H llc90 v4
%             sNx =  30,
%             sNy =  30,


% conver the XC YC coordinates to patchface.
xgrid = patchface3D(nx, nx*13, 1, XC , 2);
ygrid = patchface3D(nx, nx*13, 1, YC , 2);

list_in ={'xgrid','ygrid','XC11','YC11','XCNINJ','YCNINJ','iTile','jTile','tileNo'};
list_out={'lon'  ,'lat'  ,'XC11','YC11','XCNINJ','YCNINJ','i'    ,'j'};

tileCount=0;

%get 5 faces
[temp,xgrid]=patchface3D(4*nx,4*nx,1,xgrid,0);clear temp
[temp,ygrid]=patchface3D(4*nx,4*nx,1,ygrid,0);clear temp

% initalize some variables
XC11=xgrid; YC11=ygrid; XCNINJ=xgrid; YCNINJ=xgrid; iTile=xgrid; jTile=xgrid; tileNo=xgrid;

% loop through each of the 5 faces in the grid structure;
for iF=1:size(xgrid,2);
    
    % pull the XC and YC
    face_XC=xgrid{iF};
    face_YC=ygrid{iF};
    
    [iF]

    % loop through each tile - of which there are size(face_XC,1) /ni in i
    % and size(face_XC,2)/nj in j
    for ii=1:size(face_XC,1)/ni;
        for jj=1:size(face_XC,2)/nj;
            % accumualte tile counter
            tileCount=tileCount+1;

            % find the indicies of this particular tile in i and j
            tmp_i=[1:ni]+ni*(ii-1);
            tmp_j=[1:nj]+nj*(jj-1);
            
            % pull the XC and YC of this tile
            tmp_XC=face_XC(tmp_i,tmp_j);
            tmp_YC=face_YC(tmp_i,tmp_j);

            % pull the XC and YC at position 1,1 of this tile
            XC11{iF}(tmp_i,tmp_j)=tmp_XC(1,1);
            YC11{iF}(tmp_i,tmp_j)=tmp_YC(1,1);

            % pull the XC and YC at the position (end, end) of this tile
            XCNINJ{iF}(tmp_i,tmp_j)=tmp_XC(end,end);
            YCNINJ{iF}(tmp_i,tmp_j)=tmp_YC(end,end);
            
            % fill in this iTile/jTile structure at this tile points with a
            % local count of the i and j within the tile [from 1:ni] 
            % basically, this is a map of where in the tile each i,j point
            % is.  in llcv4 it is 1:30 in i and 1:30 in j            
            iTile{iF}(tmp_i,tmp_j)=[1:ni]'*ones(1,nj);
            jTile{iF}(tmp_i,tmp_j)=ones(ni,1)*[1:nj];

            % this is a count of the tile.
            tileNo{iF}(tmp_i,tmp_j)=tileCount*ones(ni,nj);
        end;
    end;
end;

%  now take these funky structures, cast them into patchface form, then use
%  prof point to pull the value at the profile point that we need
for k=1:size(list_out,2);
    % puts this in the original llc fucked up face
    eval(['temp=' list_in{k} ';']);temp=[temp{1},temp{2},temp{3},temp{4}',temp{5}'];

    % puts this into patchface format
    eval([list_in{k} '=patchface3D(nx,13*nx,1,temp,3);']);
    
    % use the prof_point to pull the right value from whatever list_in{k}
    % is.. list_in{k} is in patchface format, from above.    
    eval(['MITprof.prof_interp_' list_out{k} '=' list_in{k} '(MITprof.prof_point);']);
end;

%one last thing: "weights", which is 1 b/c we're using nearest neighbor:
MITprof.prof_interp_weights = ones(size(MITprof.prof_point));


