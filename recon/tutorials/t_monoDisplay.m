
% t_monoDisplay
%
% Create a display with monochromatic primaries, and demonstrate
% transformation from linear coordinates in a more standard diplsay
% to this one.

%% Clear
clear; close all;

%% Load up a standard display
%displayFile = 'CRT12BitDisplay.mat';
%display = load(fullfile(dataBaseDir, displayFile));
origDisplayFile = 'LCD-Apple';
origDisplay = displayCreate(origDisplayFile);
wls = displayGet(origDisplay,'wave');
S = WlsToS(wls);
if (S(2) > 4)
    error('Want finer wavelength spacing to represent monochromatic primaries');
end

%% Load in cone spectral sensitivities
coneData = load('T_cones_ss2');
T_cones = SplineCmf(coneData.S_cones_ss2,coneData.T_cones_ss2,S);

%% Load in an RGB image
imageSettings = im2double(imread('eagle.jpg'));
[scene,imageSettingsCheck,imagePrimary] = sceneFromFile(imageSettings, 'rgb', [], origDisplay);
if (max(abs(imageSettings(:)-imageSettingsCheck(:))) ~= 0)
    error('Scene from file not behaving as expected');
end
vcAddObject(scene);
sceneWindow;

%% The linearized image is what we call a primary image in PTB
[imagePrimaryCal,imageM,imageN] = ImageToCalFormat(imagePrimary);

%% Get the display calibration information in a form DHB understands
gammaMethod = 1;
origCalStruct = ptb.GeneratePTCalStructFromIsetbioDisplayObject(origDisplay);
origCalObj = ObjectToHandleCalOrCalStruct(origCalStruct);
SetSensorColorSpace(origCalObj,T_cones,S);
SetGammaMethod(origCalObj,gammaMethod);
if (any(origCalObj.get('P_ambient') ~= 0))
    error('Below we assume display ambient is zero.');
end

%% Get radiance using PTB knowledge
origPDevice = origCalObj.get('P_device');
origRadiancePTBCal = origPDevice*imagePrimaryCal;
origExcitationsPTBCal = T_cones*origRadiancePTBCal;

%% Get the radiance image out of the ISETBio scene, and into PTB convention.
%
% Verify this works as expected
origRadianceFromScence = sceneGet(scene,'energy')*S(2);
origRadianceFromScenceCal = ImageToCalFormat(origRadianceFromScence);
origExcitationsFromSceneCal = T_cones*origRadianceFromScenceCal;
if (max(abs(origRadiancePTBCal(:)-origRadianceFromScenceCal(:))/mean(origRadianceFromScenceCal(:))) > 1e-6)
    error('Cannot connect PTB and ISETBio worlds the way we should be able to');
end



% %%
% ISETBioDisplayObject = ptb.GenerateIsetbioDisplayObjectFromPTBCalStruct('SACC', screenCalStruct, extraCalData, false);
% ISETBioDisplayObject = rmfield(ISETBioDisplayObject,'dixel');

