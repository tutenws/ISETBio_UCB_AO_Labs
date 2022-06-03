function [xValsVector, yValsVector, rfSpacingDegs] = generateHexMosaic(theConeMosaic, filterSizeMultiplier)

% This function will generate a haxagonal grid with element spacing defined
% relative the cone spacing in an ISETBio cMosaic object.
%
% INPUTS
%  theConeMosaic            cMosaic object from ISETBio.
%  filterSizeMultiplier     The variable defines the spacing between
%                           elements in the hex array as a multiple of the
%                           mean cone spacing in the cMosaic object.
%                           Default value is 3.
%
%
% OUTPUTS
%  xValsVector              x-coordinates of the hex array, in retinal
%                           coordinates (degrees) defined by the cMosaic
%                           object.
%  yValsVector              y-coordinates as above.
%
% HISTORY
%  06-02-2022               wst wrote it.% 

if nargin < 2
    filterSizeMultiplier = 3;
end
if nargin < 1
    error('No cone mosaic object specified.');
end

% Based on this number, what is the col (i.e. x) and row (i.e. y) spacing
% of the hex array?
colSpacing = filterSizeMultiplier.*mean(theConeMosaic.Mosaic.coneRFspacingsDegs);
rfSpacingDegs = colSpacing;
rowSpacing = (sqrt(3)/2).*colSpacing;

% Get the mosaic center from theConeMosaic object
yCenter = theConeMosaic.Mosaic.eccentricityDegs(2);
xCenter = theConeMosaic.Mosaic.eccentricityDegs(1);

% Compute the number of rows and colums, with a bit of padding; consider
% how to handle post-receptoral summation units near the edge of the mosaic
% sometime later
numRows = round(theConeMosaic.Mosaic.sizeDegs(2)./rowSpacing) + 2;
numCols = round(theConeMosaic.Mosaic.sizeDegs(1)./colSpacing) + 2;

% Initialize where the rows should stop and start; the initial grid will be
% shifted to align with the center of the cone mosaic later.
rowStart = yCenter - theConeMosaic.Mosaic.sizeDegs(2)./2 - rowSpacing;
colStart = xCenter - theConeMosaic.Mosaic.sizeDegs(1)./2 - colSpacing;

% Generate a grid of X and Y values
xValVectorSeed = colStart:colSpacing:(colSpacing*numCols)+colStart;
yValVectorSeed = rowStart:rowSpacing:(rowSpacing*numRows)+rowStart;
[xVals, yVals] = meshgrid(xValVectorSeed, yValVectorSeed);

% Shift every other row by half the column spacing so that the grid becomes
% hexagonal
xVals(2:2:end,:) = xVals(2:2:end,:)+colSpacing/2; 

% Reformat into a vector
xValsVector = reshape(xVals, numel(xVals), 1);
yValsVector = reshape(yVals, numel(yVals), 1);

% Figure out the x- and y- coordinates nearest the mosaic center
[centerIDX, ~] = knnsearch([xValsVector yValsVector],[xCenter yCenter]);

% Compute the shift so the hex grid aligns with the mosaic center
xShift = xValsVector(centerIDX) - xCenter;
yShift = yValsVector(centerIDX) - yCenter;

% Apply shift
xValsVector = xValsVector-xShift;
yValsVector = yValsVector-yShift;