% Housekeeping
clear all;
close all;
clc;

% Load in a data file
load('C:\Users\William Tuten\Documents\MATLAB\projects\ISETBio_UCB_AO_Labs\coneGainModel\results\propL_0p47_2_arcmin.mat');

% Show the cone mosaic object
theConeMosaic.Mosaic.visualize;
hold on;

% Generate the post-receptoral array
postReceptoralSpacingConeSpacingMultiplier = 3;
[xPRA, yPRA, rfSpacingDegs] = generateHexMosaic(theConeMosaic, postReceptoralSpacingConeSpacingMultiplier);

% For each cone in the mosaic, find the distance (coneDistToPRA_degs) to
% the nearest element in the post-receptoral array (idxPRA)
[idxPRA, coneDistToPRA_degs] = knnsearch([xPRA yPRA], theConeMosaic.Mosaic.coneRFpositionsDegs);

% For each post-receptoral element, determine how many cones are in its
% "receptive field"
coneCountsPerPRA = zeros(size(xPRA));
for n = 1:length(xPRA)
    coneCountsPerPRA(n) = length(idxPRA(idxPRA==n));
end

% For each post-receptoral element, determine the underlying cone ratio
proportionLPerPRA = nan(size(xPRA)); 
rfType = 'circular';
for n = 1:length(xPRA)
    coneIndices = find(idxPRA==n);
    coneTypes = theConeMosaic.Mosaic.coneTypes(coneIndices);
    coneDistances = coneDistToPRA_degs(coneIndices);
    if ~isempty(coneIndices)
        switch rfType
            case 'gaussian'
                % For now, let's make the FWHM equal to the element spacing
                gaussianFWHM = rfSpacingDegs;
                gaussianSTD = gaussianFWHM./2.35482;
                lConeWeights = exp(-((coneDistances(coneTypes==1).^2)/(2*(gaussianSTD.^2))));
                mConeWeights = exp(-((coneDistances(coneTypes==2).^2)/(2*(gaussianSTD.^2))));
                proportionLPerPRA(n) = sum(lConeWeights)./(sum(lConeWeights) + sum(mConeWeights));
            case 'circular' % No distance-dependent weighting
                proportionLPerPRA(n) = sum(coneTypes==1)./length(coneTypes(coneTypes~=3));
            otherwise 
                error('Improper RF type specified.')
        end
    else
        % Do nothing
    end
end