% First stab at modeling some unique yellow experiments in ISETBio
%
% The goal of this exercise is to develop a framework for modeling unique
% yellow experiments in ISETBio. We'll build from routines outlined in some
% of the ISETBio livescripts. To start, we will model the L- and M-cone
% isomerizations in a phsyiologically-plausible cone mosaic. The stimuli we
% use will be monochromatic wavelengths centered around the conventional
% unique yellow setting (578 nm). Once we know the L and M isomerizations
% caused by 578 nm light, we'll simulate a 2 AFC experiment where the
% observer's task is to report whether the light is reddish or greenish
% based on the 578-nm criterion. Eventually, we'll want to run this for a
% two-primary experiment, with and without AO correction, and with and
% without compensation for fixational eye movements, and for various
% stimulus sizes and durations.
%
% We also need to sort out how to build up the optics model so that we can
% simulate AO-corrected stimuli that have been adjusted to correct for the
% eye's LCA. At present I do not know how to accomplish this.
%
% 6/27/20       wst getting started on this

%% Housekeeping
clear all
close all
clc;

%% Step 0: Initial stimulus parameters
% ISETBio represents stimuli as "scenes". Most of the plug-and-play
% routines involve monitor-based stimuli, so we have to hack something
% together as best we know how.

% What wavelengths are at play here
stimParams.wls = [543 578 680]; % Monochromatic wavelengths
stimParams.sceneSizePixels = 256; % Number of pixels in the raster; may want to downsample this to save computation time
stimParams.sceneSizeDeg = 0.5; % Size of the scene, in degrees; should be at least as large as the largest stimulus
stimParams.scenePPD = stimParams.sceneSizePixels/stimParams.sceneSizeDeg; % Define scaling here

% Make a circular mask to transform the initial scene from uniform square
% to a stimulus
stimParams.spotDiamArcmin = 7.5;
stimParams.spotDiamDeg = stimParams.spotDiamArcmin./60;
stimParams.spotDiamPixels = stimParams.spotDiamDeg.*stimParams.scenePPD;
stimParams.spotMask = double(Circle(stimParams.spotDiamPixels/2)); % Make the mask using PTB Circle
% Pad mask to match scene size
borderSizePixels = 0.5.*(stimParams.sceneSizePixels-size(stimParams.spotMask,1));
stimParams.spotMask = padarray(stimParams.spotMask, [floor(borderSizePixels) floor(borderSizePixels)], 0, 'pre');
stimParams.spotMask = padarray(stimParams.spotMask, [ceil(borderSizePixels) ceil(borderSizePixels)], 0, 'post');

if size(stimParams.spotMask,1) ~= stimParams.sceneSizePixels
    error('Scene and mask sizes are unequal');
end

%% Step 1: Create the cone mosaic
% Generate a hexagonal cone mosaic with ecc-based cone quantal efficiency
propS = 0.1; % Proportion S-cones (10% is default for ISETBio, so we'll use that here)
LM_ratio = 1; % Ratio of L to M cones
propL = (1-propS).*(LM_ratio./(LM_ratio+1)); % Convert to proportion L
propM = 1-propS-propL;  % Proportion M
coneMosaicDist = [0 propL propM propS]; % This gets input into the below routine
theMosaic = coneMosaicHex(5, ...              % hex lattice sampling factor
    'fovDegs', stimParams.sceneSizeDeg, ...        % match mosaic width to stimulus size
    'eccBasedConeDensity', true, ...           % cone density varies with eccentricity
    'eccBasedConeQuantalEfficiency', true, ... % cone quantal efficiency varies with eccentricity
    'integrationTime', 30/1000, ...            % 30 msec integration time (1 AOSLO frame)
    'spatialDensity', coneMosaicDist, ...      % [something propL propM propS]
    'maxGridAdjustmentIterations', 50);        % terminate iterative lattice adjustment after 50 iterations


%% Step 2: Generate the scene at each wavelength

% Intial creation of the scene, which will produce equal energies at each
% wavelength in the "wls" variable; the size of the scene is sceneSize x
% sceneSize x length(wls); ADD LOOP
for w = 1:length(stimParams.wls)
    sceneSpot = sceneCreate('uniform monochromatic', stimParams.wls(w), stimParams.sceneSizePixels);
    
    % Set the visual angle and viewing distance of the scene
    sceneSpot = sceneSet(sceneSpot, 'wAngular', stimParams.sceneSizeDeg, 'distance', 0.57);
    
    % Adjust the luminance of each wavelength layer to be the same
    

    
    % Multiply photon matrix in scene by stimulus mask; then update scene
    sceneSpot = sceneSet(sceneSpot, 'photons', sceneSpot.data.photons.*stimParams.spotMask);
    
    
    %% Step 3: Create the diffraction-limited optical image with no LCA
    pupilDiameterMm = 6;
    accommodatedWavelength = stimParams.wls(1); % Pin best focus here (usually will be 543 nm)
    zCoeffs = zeros(66,1); % Simulating diffraction-limited performance here (add flag later to compare to ordinary aberrations)
    
    % Create wvfP object
    wvfP = wvfCreate('calc wavelengths', stimParams.wls, 'zcoeffs', zCoeffs, ...
        'name', sprintf('human-%d', pupilDiameterMm));
    wvfP = wvfSet(wvfP, 'measured pupil size', pupilDiameterMm);
    wvfP = wvfSet(wvfP, 'calc pupil size', pupilDiameterMm);
    % Deal with best focus by specifying that the wavefront parameters
    % were measured at the wavelength we want to say is in focus. This
    % is a little bit of a hack but seems OK for the diffraction limited case
    % we're using here.
    wvfP = wvfSet(wvfP, 'measured wavelength', accommodatedWavelength);
    % Compute pupil function using 'no lca' key/value pair to turn off LCA.
    % You can turn it back on (i.e. to "false") to compare the effect.
    wvfPNoLca = wvfComputePupilFunction(wvfP,false,'no lca',true); % First "false" is for "showBar"?
    wvfPNoLca = wvfComputePSF(wvfPNoLca); % Update wvfPNoLCA object with PSF
    theOI = wvf2oi(wvfPNoLca); % Covert to wvfP to OI
    opticsNoLca = oiGet(theOI, 'optics'); % Pull in the optics
    opticsNoLca = opticsSet(opticsNoLca, 'model', 'shift invariant'); % Update the optics model
    opticsNoLca = opticsSet(opticsNoLca, 'name', 'human-wvf-nolca'); % Rename for clarity
    theOI = oiSet(theOI,'optics',opticsNoLca, 'wAngular', stimParams.sceneSizeDeg); % Update the OI
    
    % Visualize OI
    visualizedSpatialSupportArcMin = 6; % Size of displayed PSF, in arc min
    visualizedSpatialSfrequencyCPD = 120; % Max spatial frequency for MTF
    % Visualize the PSF/OTF at the shortest wavelength (diffraction-limited)
    visualizeOptics(theOI, stimParams.wls(1), visualizedSpatialSupportArcMin, visualizedSpatialSfrequencyCPD);
    % Visualize the PSF/OTF at the longest wavelength (also diffraction-limited, but MTF
    % should be a little worse "because physics")
    visualizeOptics(theOI, stimParams.wls(2), visualizedSpatialSupportArcMin, visualizedSpatialSfrequencyCPD);
    
    
    %% Step 4: Compute the retinal image
    % First, compute the retinal image
    theOI = oiCompute(theOI, sceneSpot);
    % Visualize different aspects of the computed optical image
    visualizeOpticalImage(theOI, 'displayRetinalContrastProfiles', true);
    
    %% Step 5: Set up eye movement paths
    % Generate 2 instances of zero movement eye movement paths
    nTrialsNum = 10; % This times the integration time of "theMosaic" gives stimulus duration? Need to learn more
    emPath = zeros(nTrialsNum, 1, 2); % Zeros for now; will want to incorporate later
    
    %% Step 6: Compute mosaic isomerizations
    coneExcitations = theMosaic.compute(theOI, 'emPath', emPath);
    % Visualize cone mosaic and its cone excitation responses
    visualizeConeMosaicResponses(theMosaic, coneExcitations, 'R*/cone/tau');
end