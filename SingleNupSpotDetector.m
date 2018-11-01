%% SingleNupSpotDetector

% Ideally, this code will, upon being given an image and the proposed
% coordinates of a single Nup spot, determine whether or not the spot is
% valid (by applying a couple of filtering tests) and then compute the
% integrated intensity of the spot by applying a Gaussian fit. Filtering
% tests will be size and eccentricity. Additionally, have the code check
% the brightness of the proposed spot in the planes above and below the
% plane in which I identified it. If the spot is brighter in one of those
% two planes, then that suggests that the plane I was looking at isn't the
% focal plane for that spot.

clc
clear all

% Images used with this program must be stored in the MATLAB home folder.
% This is easily modifiable if desired.

%% Configuration Variables
% The plus one is there to adjust for the fact that arrays in
% MATLAB don't start at the same values as FIJI arrays! I think it's pretty
% important to experiment with the crop radius and ways to compute an
% appropriate crop radius, because it does seem to affect the peak
% intensity quite a bit. I'd like to programmatically determine both the
% cropRadius and the backgroundIntensity.

fileID = fopen([pwd, '/config.txt'],'r');
imageNamePrecursor = fgetl(fileID);
spacesLocatedAt = find(imageNamePrecursor == ' ');
imageName = imageNamePrecursor(spacesLocatedAt(2) + 1 : size(imageNamePrecursor, 2));

spotCoordinatesPrecursor = string(fgetl(fileID));
spotCoordinates = split(spotCoordinatesPrecursor, ",");
planeNumber = double(spotCoordinates(1,1));
x_coord = double(spotCoordinates(2,1)) + 1;
y_coord = double(spotCoordinates(3,1)) + 1;

cropRadius = 4;
cropDiameter = 2*cropRadius + 1;
backgroundIntensity = 400;
spotID = 0;

%% Image Read and Spot Identification
I = imread([pwd , '/', imageName , '.tif'], planeNumber);

% Scales a copy of the image in order to provide a visual reference. Should
% only be used during testing. Create a parameter (from config file) that enables/disables
% this section of code?
I_scaled = imadjust(I);
figure(1)
imshow(I_scaled, 'InitialMagnification', 'fit');
hold on
plot(x_coord, y_coord, 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 20);

%% Crop Image Around Identified Spot
% Scales a copy of the cropped image (an 9x9 box around the proposed
% spot) in order to provide a visual reference. Use local adaptive
% thresholding to identify the spot, and then crop the image to the spot's
% bounding box. 
croppedImage = I((y_coord - cropRadius) : (y_coord + cropRadius), (x_coord - cropRadius) : (x_coord + cropRadius));
doubleCroppedImage = cast(croppedImage, 'double');
figure(2);
imshow(imadjust(croppedImage), 'InitialMagnification', 'fit');

%% Fit Cropped Image to 2D Gaussian Curve
xdata = zeros(2, cropDiameter.^2);

for i=1:cropDiameter ^ 2
    if(rem(i, cropDiameter) ~= 0)
        xdata(1, i) = rem(i, cropDiameter);
    else
        xdata(1, i) = cropDiameter;
    end
end

for i=1:cropDiameter
    xdata(2, (cropDiameter * i - (cropDiameter - 1)) : (cropDiameter * i)) = i;
end

% I'm still a bit sketched out here and concerned that the dimensions might
% be going in the wrong order. Hopefully that isn't the case.
ydata = reshape(doubleCroppedImage, [1, cropDiameter.^2]);

% predicted is an anonymous fitting function that lsqcurvefit will fit the
% data to. a will be a vector with five elements: the amplitude, the x shift, the x
% standard deviation, the y shift, and the y standard deviation. 

predicted = @(a, xdata) a(1) * exp(-((xdata(1, :) - a(2)).^2 / (2 * (a(3).^2))) - ((xdata(2, :) - a(4)).^2) / (2 * (a(5).^2)));

% a0 is the first estimate of parameters that the lsqcurvefit will use. 
a0 = [doubleCroppedImage(cropRadius+1,cropRadius+1); 0; 4; 0; 4];

% Performs the curve fitting.
opts = optimset('Display','off');
[ahat, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(predicted, a0, xdata, ydata, [], [], opts);

% The final fitted Gaussian function. Can be used for the integration.
gaussianFit = reshape(predicted(ahat, xdata), [9,9]);
thresholdedFit = imbinarize(gaussianFit, backgroundIntensity);

%% Compute Properties of the 2D Gaussian Fit
% Use regionprops to extract the size of the portion that is above the
% threshold?
eccentricity = sqrt(1 - min( (ahat(3).^ 2 / ahat(5).^ 2), (ahat(5).^ 2 / ahat(3).^ 2) ) );

% Compute the integrated intensity over a given set of bounds (which really should be computed programatically):
% xbounds = [-1, 1];
% ybounds = [-1, 1];

% I'm skeptical of this for the reason that I should really be integrating
% over an ellipse rather than over a square.
% intIntensity = ahat(1) * (pi/2) * ahat(3) * ahat(5) * ( erf(xbounds(2) / (ahat(3) * sqrt(2))) - erf(xbounds(1) / (ahat(3) * sqrt(2)))) * ( erf(ybounds(2) / (ahat(5) * sqrt(2))) - erf(ybounds(1) / (ahat(5) * sqrt(2))));

% This more primitive method simply sums up all of the pixel values that
% exceed the background. Then it subtracts the backgroundIntensity,
% multiplied by the number of pixels that exceed the background
unintIntensity = sum(sum(thresholdedFit .* gaussianFit)) - backgroundIntensity * sum(sum(thresholdedFit));

% plots the centre of the proposed Gaussian (which should be the centre of
% the cropped image). NTS: double-check that this is working right!
hold on
plot(ahat(4), ahat(2), 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 20);

figure(3);
heatmap(gaussianFit);

%% Output Data into a .csv File
% Note that files being opened for writing should be opened with 'w'
% permission, which will delete the previous contents, or with 'a'
% permission, which will append new text to the previous contents.
fileID = fopen([pwd, '/SpotDetectorOutput.csv'],'w');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'spotID', 'plane', 'x', 'y', 'peak', 'x0', 'xdev', 'y0', 'ydev', 'ecc', 'int');
fprintf(fileID, '%u,%u,%u,%u,%5.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.3f,%5.2f', spotID, planeNumber, x_coord, y_coord, ahat(1), ahat(2), ahat(3), ahat(4), ahat(5), eccentricity, unintIntensity);
fclose(fileID);

