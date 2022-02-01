% Make AO scene for unique yellow experiments
%
% Description:
%    First-pass attempt at setting up a ISETBio scene that
%    approximates the conditions for the small-spot, retinally stabilized
%    unique yellow experiments by Boehm et al.
% 
% See also: t_AORectangleScene

% History: 
%   01/28/22  wst  Wrote it form t_AORectangleScene.
%   01/29/22  dhb  Tune up.

%% Clear and close
clear; close all;

%% Set parameters 
defocusDiopters = 0.05;
pupilDiameterMm = 7;
pupilAreaMm2 = pi*(pupilDiameterMm/2)^2;
 
% Rectangle horizontal and vertical size in arcminutes. (WST updated)
spotHeightMinutes = 2.2;
spotWidthMinutes= 2.2;

% Position of rectangle. These numbers give the center position in arcmin.
%    Vertical - positive is up 
%    Horizontal - positive is right
spotVerticalPositionMinutes = 2.5;
spotHorizontalPositionMinutes = 5;

% Other stimulus parameters
%
% Define basic parameters of the AO stimulus
wls = (400:10:750)';
S = WlsToS(wls);

% For now, let's keep our conditions monochromatic rather than trying to
% model what we actually did, which was to make yellow with a mixture of
% our red (680 nm) and green (543 nm) primaries
spotWavelengthNm = 580; 
spotFWHMNm = 5;

% Stimulus luminance, according to Ally, was around 2844 cd/m2. The
% following lines try to work from that number back out to a sensible
% corneal irradiance, in uW, for 580 nm light
spotLuminanceCdM2 = 2844;

% Spatial parameters
nPixels = 128;
fieldSizeMinutes = 60;
fieldSizeDegs = fieldSizeMinutes/60;
fieldSizeDegs2 = fieldSizeDegs^2;



% Create a relative spectrum where the only nonzero value is at the spot
% wavelength
relativeSpectrumSpot = zeros(size(wls));
wlInd = find(wls==spotWavelengthNm);
relativeSpectrumSpot(wlInd) = 1;  

% Go from lumninance to radiance. Radiance comes back using PTB convention
% of power per wlband, not power per nm.
[spotRadianceWattsM2SrWlband, spotRadianceS] = LumToRadiance(relativeSpectrumSpot, S, spotLuminanceCdM2, 'Photopic');

% And then go from radiance back out to corneal irradiance given our
% field size
cornealIrradianceWattsM2Spot = RadianceAndDegrees2ToCornIrradiance(spotRadianceWattsM2SrWlband, fieldSizeDegs2);

% Pass the result at 580 nm in uW to the scene generation code below
Mm2ToM2 = 1e-6;
WToUW = 1e6;
spotIrradianceUWM2Spot = cornealIrradianceWattsM2Spot(wlInd) * WToUW;
spotPowerUW = spotIrradianceUWM2Spot*pupilAreaMm2*Mm2ToM2;

% Let's convert spotPower to nW/deg^2, about which DHB has some intutition.
% That is, we use about 70 nW/deg^2 in the Penn AOSLO for full power
% visible stimulus at 545 nM or so.  The number here comes out ballbark
% reasonable.
UWToNW = 1000;
spotPowerNwDeg2 = spotPowerUW*UWToNW/fieldSizeDegs2;

% The underlying routine starts with a background and works in contrast.
% To deal with this, we specify a background down from desired power by
% some large factor and then take contrast as this factor.  This produces a
% spot with desired power and a background of essentially 0.  We'll patch
% that up later - doing it this way interfaces to the extant routine and is
% a kluge, but one I think we can live with.
bgFactor = 5000;
spotBgPowerUW = spotPowerUW/bgFactor;

% Convert our parameters to those that matter to the two spot routine.
% Using stimAngle of 90 means we only see spot 2, which is drawn second and
% thus always over spot 1.
stimAngle = 90;

% This is center-to-center separation of the two spots. So with values of
% zero they overlap. Because the separation is accomplished by applying
% half of it to each underlying spot, we multiply by two to convert
% position offsets to separations.  The minus sign makes the sign
% convention as described above.
spotVerticalSepMinutes = -2*spotVerticalPositionMinutes;
spotHorizontalSepMinutes = -2*spotHorizontalPositionMinutes;

% Set up two spot params
rectParams = struct(...
    'type', 'basic', ...                            % type
    'viewingDistanceMeters', 2, ...                 % viewing distance: in meters
    'wls', 400:10:750, ...                          % wavelength support of the primary SPDs: in nanometers
    'stimAngle', stimAngle, ...                     % stimulus angle in incr/decr plane
    'spotWl', spotWavelengthNm, ...                 % spot wavelength: in nm
    'spotFWHM', spotFWHMNm, ...                     % spot full width at half max: in nm
    'spotWidthDegs', spotWidthMinutes/60, ...       % spot width: in degrees
    'spotHeightDegs', spotHeightMinutes/60, ...     % spot height: in degrees
    'spotVerticalSepDegs', spotVerticalSepMinutes/60 , ... % spot center to center vertical separation: in degrees
    'spotHorizontalSepDegs', spotHorizontalSepMinutes/60, ... % spot center to center horizontal separation: in degrees
    'spotBgDegs', fieldSizeDegs, ...                % spot background: in degrees
    'spotBgPowerUW', spotBgPowerUW, ...             % spot background power: in uW
    'imagingWl', 750, ...                           % imaging beam wavelength: in nm
    'imagingFWHM', 20, ...                          % imaging beam full width at half max: in nm
    'imagingBgPowerUW', 0, ...                      % imaging beam power entering eye: in uW
    'fovDegs', fieldSizeDegs, ...                   % spatial: full scene field of view, in degrees
    'pixelsNum', nPixels, ...                       % spatial: desired size of stimulus in pixels
    'temporalModulation', 'flashed', ...            % temporal modulation mode: choose between {'flashed'}
    'temporalModulationParams', struct(...          % temporal: modulation params struct
      'stimOnFrameIndices', [1], ...                %   params relevant to the temporalModulationMode
      'stimDurationFramesNum', 1), ...              %   params relevant to the temporalModulationMode
    'frameDurationSeconds', 3/16, ...               % temporal: frame duration, in seconds
    'pupilDiameterMm', pupilDiameterMm ...          % pupil diameter mm
    );
                     
% Create a static two-spot AO scene with a particular incr-decr direction,
% and other relevant parameters. This uses compute function handle for
% two-spot AO stimuli, as commented above the parameters have been cooked
% to make one rectangular spot on a dark background.
rectComputeFunction = @sceAOTwoSpot;

% Instantiate a sceneEngine with the above sceneComputeFunctionHandle
% and the custom params.
rectSceneEngine = sceneEngine(rectComputeFunction, rectParams);

% Create a scene sequence at a particular contrast.  Here there is just one
% frame in the returned cell array of scenes. We pull out that first entry.
spotContrast = bgFactor*1.0;
rectSceneSequence = rectSceneEngine.compute(spotContrast);
rectScene = rectSceneSequence{1};

% Adjust the background.  First we'll zero out the very
% small numbers we put in to be able to specify spot power in the form of a
% contrast.
sceneRadiancePhotonsM2SrSecNm = sceneGet(rectScene,'photons');
maxPhotons = max(sceneRadiancePhotonsM2SrSecNm(:));
sceneRadiancePhotonsM2SrSecNm(sceneRadiancePhotonsM2SrSecNm(:) < 4*maxPhotons/bgFactor) = 0;

% Now approximate the background, which is effectively metameric to equal
% energy white and has a luminance of 26 cd/m2
bgLuminanceCdM2 = 26;

% Relative spectrum reflects equal energy at all wavelengths (not true, but ok for now)
relativeSpectrumBg = ones(size(wls)); 
bgRadianceWattsM2SrSecWlband = LumToRadiance(relativeSpectrumBg, S, bgLuminanceCdM2, 'Photopic');

% Convert to photonsM2SrSecNm, which are the units that ISETBio expects
bgRadiancePhotonsM2SrSecWlband = EnergyToQuanta(wls,bgRadianceWattsM2SrSecWlband);
bgRadiancePhotonsM2SrSecNm = bgRadiancePhotonsM2SrSecWlband/S(2);

% Add background to spot image
for ww = 1:length(wls)
    sceneRadiancePhotonsM2SrSecNm(:,:,ww) = ...
        sceneRadiancePhotonsM2SrSecNm(:,:,ww) + bgRadiancePhotonsM2SrSecNm(ww);
end

% Push adjusted image back into the scene. Note that default ISETBio scene
% is not represented at full double precision.  Setting bit depth to 64
% uses double precision which is what we want here.
rectScene = sceneSet(rectScene, 'bit depth', 64);
rectScene = sceneSet(rectScene,'photons',sceneRadiancePhotonsM2SrSecNm);
sceneRadianceCheck = sceneGet(rectScene,'photons');
if (max(abs(sceneRadianceCheck(:) - sceneRadiancePhotonsM2SrSecNm(:)))/mean(sceneRadiancePhotonsM2SrSecNm(:)) > 1e-10)
    error('Failure in scene set');
end

% Visualize.
%
% Use the Analyze:Line (luminance):/Horizontal menu item to select a line
% and see plot of luminance profile through the center.  Both spot and bg
% luminace come out as specified (remember that spot luminance is the sum
% of bg and specified spot luminances).
%
% Use the Plot:Radiance (photons) menu to select and ROI and plot mean
% spectrum in that region. Selecting background shows an equal energy
% spectrum, while selecting spot shows monochromatic spot on a relatively
% weak equal energy background. Not that an equal energy spectrum is not an
% equal quantal spectrum, though. You can verify equal energy by using the
% option that shows (energy) rather than photons.
ieAddObject(rectScene);
sceneWindow;