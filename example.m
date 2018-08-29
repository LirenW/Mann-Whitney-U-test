%%
% This progrom is for the Mann Whitney U (MW-U) test. The main 
% function in MWtest.m. Further descriptions on MWtest.m
%
%  Created by Liren
%  May 2018
%
%%

clc
clear
tic
    N       = 1000;                     % Pairs of bootstrap resampling
    flag    = false;                    % flag - Use bootstrap or not. True for open the bootstrap
                                        % resampling process.
    fieldA  = './';  % Two input field
    fieldB  = './';
    outFile = './example.nc';
    varName = 'var';

    %% Reading process

    lon     = ncread(fieldA, 'lon');
    lat     = ncread(fieldB, 'lat');
  
    varA    = ncread(fieldA, varName);
    varB    = ncread(fieldB, varName);

    dimSize = size(varA);
    NLON    = dimSize(1);
    NLAT    = dimSize(2);
    NTIME   = dimSize(3);
    
    %% Calculation process

    P       = MWtest(varA ,varB, flag, N, NLAT, NLON, NTIME);

    %% Write the netCDF file

    ncid     = netcdf.create(outFile , 'NETCDF4');

    % Global Attributions setting
    varidGlo = netcdf.getConstant('GLOBAL');
    netcdf.putAtt(ncid , varidGlo , 'creation_date' , datestr(now));
    if flag == true
        netcdf.putAtt(ncid , varidGlo , 'bootstrap_sample' , N);
    else
        netcdf.putAtt(ncid , varidGlo , 'bootstrap_sample' , 0);
    end

    % Latitude & Longitude Attributions setting
    dimidx   = netcdf.defDim(ncid , 'lat' , NLAT);
    dimidy   = netcdf.defDim(ncid , 'lon' , NLON);

    varidLat = netcdf.defVar(ncid , 'lat' , 'float' , dimidx);
    varidLon = netcdf.defVar(ncid , 'lon' , 'float' , dimidy);

    netcdf.putAtt(ncid , varidLat , 'units' , 'degrees_north');
    netcdf.putAtt(ncid , varidLon , 'units' , 'degrees_east' );

    varidVar = netcdf.defVar(ncid , 'P' , 'double' , [dimidy dimidx] );
    netcdf.putAtt(ncid , varidVar , 'long_name' , 'p-value' );

    netcdf.defVarFill(ncid , varidVar , false , -999); % Must under netCDF4

    netcdf.endDef(ncid);

    netcdf.putVar(ncid , varidLat  , lat);
    netcdf.putVar(ncid , varidLon  , lon);
    netcdf.putVar(ncid , varidVar  , P  );

    netcdf.close(ncid);

toc