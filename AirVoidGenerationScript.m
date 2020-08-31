%% AirVoidGenerationScript uses the GenerateCircleLogNormalDistribution function to generate circular voids.
%
%  Discussion:
%
%    The GenerateCircleLogNormalDistribution code is a function 
%    that takes the Matrix M as an input, so first we must determine 
%    the Matrix M which contains the coordinates of all aggregate 
%    particles ordered with respect to their coordinates.
%
%  Copyright:
%
%    If you want to use these subroutines (research ONLY), please cite our papers:
%   
%  Impact of Void Morphology on the Mechanical Response of Time-Dependent Heterogeneous Media: A Numerical Investigation Approach
%  Journal of Materials in Civil Engineering
%  2020-07 | journal-article
%  DOI: 10.1061/(ASCE)MT.1943-5533.0003252
%
%  Modified:
%
%    No recent modifications.
%
%  Author:
%
%    Aimane Najmeddine
%
  clear;
  clf;
  warning('off');

% Read the coordinates of all particle vertices (i.e., aggregates) 
% from the AUTOCAD file containing the X-Ray image of the material.
A = xlsread('COORDSALLAGG.xls');
% 

%PlotExtractedDataPoints(A);

%store the coordinates connecting each particle in order in a Matrix M;
M = FindPolygons(A);
%

%Generate the polygons and output two results: number of voids + void
%content
[counter_of_air_voids, air_void_content] = GenerateCircleLogNormalDistribution(M,100,0.1);
%
