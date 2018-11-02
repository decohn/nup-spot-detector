%% SingleNupSpotDetector

% Next up: improve the way that I calculate the background intensity.
% Following that, incorporate a loop so that the program can run on many
% proposed NPC spots at once.

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
% MATLAB don't start at the same values as FIJI arrays! Programmatically
% compute the background intensity!

fileID = fopen([pwd, '/config.txt'],'r');
imageNamePrecursor = fgetl(fileID);
spacesLocatedAt = find(imageNamePrecursor == ' ');
imageName = imageNamePrecursor(spacesLocatedAt(2) + 1 : size(imageNamePrecursor, 2));

spotCoordinatesPrecursor = string(fgetl(fileID));
spotCoordinates = split(spotCoordinatesPrecursor, ",");
planeNumber = double(spotCoordinates(1,1));
x_coord = double(spotCoordinates(2,1)) + 1;
y_coord = double(spotCoordinates(3,1)) + 1;

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

%% Crop Image Around Identified Point
% Scales a copy of the cropped image (an 19x19 box around the proposed
% spot) in order to provide a visual reference. Use local adaptive
% thresholding to identify the spot itself, and then crop the image to the spot's
% bounding box. 

% This below value for the sensitivity of the adaptive thresholding will
% really need to be played with quite a bit in order to arrive at something
% that I'm happy with. This code will binarize the entire image into a
% large number of discrete objects. I'm then going to select the single
% object that contains the spot of interest. My largest concern is how to
% optimize this thresholding. A lot of nuclei aren't nicely divided up into
% distinct NPCs, and the risk of having the bounding box be way too big
% feels reasonably large. adaptthresh also has a lot of extra optional
% parameters that I can mess with to try and improve my results.

adaptivelevel = adaptthresh(I, 0.35);
BW = imbinarize(I, adaptivelevel);
figure(2)
imshow(BW, 'InitialMagnification', 'fit');
hold on
plot(x_coord, y_coord, 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 20);

BWOI = bwselect(BW, x_coord, y_coord, 4);

% Gives the coordinates of the upper-left corner of the bounding box,
% followed by the x-width and the y-width.
object_data = regionprops(BWOI, I, {'BoundingBox'});
data_array = table2array(struct2table(object_data));
data_array(1) = cast(data_array(1), 'uint16');
data_array(2) = cast(data_array(2), 'uint16');

croppedImage = I((data_array(2)) : (data_array(2) + data_array(4) - 1), (data_array(1)) : (data_array(1) + data_array(3) - 1));

figure(3)
imshow(imadjust(croppedImage), 'InitialMagnification', 'fit');
hold on
plot(x_coord, y_coord, 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 20);

doubleCroppedImage = cast(croppedImage, 'double');

% computes the background intensity as the average of the four corners just
% outside of the bounding box. probably respectably accurate, but could
% certainly be improved further.
% fourCorners = [I(data_array(2) - 1, data_array(1) - 1), I(data_array(2) - 1, data_array(1) + data_array(3)), I(data_array(2) + data_array(4), data_array(1) - 1), I(data_array(2) + data_array(4), data_array(3) + data_array(1))];
% backgroundIntensity = mean(fourCorners);

% Instead, I chose to compute intensity by taking the average intensity of
% every pixel that is inside the bounding box but NOT inside the spot. It
% may later turn out to be better to use something like the mean of this
% value and the minimum intensity inside the bounding box. 

binaryCroppedImage = BWOI((data_array(2)) : (data_array(2) + data_array(4) - 1), (data_array(1)) : (data_array(1) + data_array(3) - 1));
backgroundRegion = cast(~binaryCroppedImage, 'uint16');
backgroundIntensity = sum(sum(backgroundRegion .* croppedImage)) / sum(sum(backgroundRegion));

%% Fit Cropped Image to 2D Gaussian Curve
xdata = zeros(2, data_array(3) * data_array(4));

for i=1:(data_array(3)*data_array(4))
    if(rem(i, data_array(4)) ~= 0)
        xdata(1, i) = rem(i, data_array(4));
    else
        xdata(1, i) = data_array(4);
    end
end

for i=1:data_array(3)
    xdata(2, (data_array(4) * i - (data_array(4) - 1)) : (data_array(4) * i)) = i;
end

% I'm still a bit sketched out here and concerned that the dimensions might
% be going in the wrong order. Hopefully that isn't the case. I think they
% are in the wrong order... there's something sketchy going on with the
% fit.
ydata = reshape(doubleCroppedImage, [1, data_array(3) * data_array(4)]);

% predicted is an anonymous fitting function that lsqcurvefit will fit the
% data to. a will be a vector with five elements: the amplitude, the x shift, the x
% standard deviation, the y shift, and the y standard deviation. 

predicted = @(a, xdata) a(1) * exp(-((xdata(1, :) - a(2)).^2 / (2 * (a(3).^2))) - ((xdata(2, :) - a(4)).^2) / (2 * (a(5).^2)));

% a0 is the first estimate of parameters that the lsqcurvefit will use. 
a0 = [double(I(x_coord, y_coord)); 0; 4; 0; 4];

% Performs the curve fitting.
opts = optimset('Display','off');
[ahat, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(predicted, a0, xdata, ydata, [], [], opts);

% The final fitted Gaussian function. Can be used for the integration.
gaussianFit = reshape(predicted(ahat, xdata), [data_array(4), data_array(3)]);
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

figure(4);
heatmap(gaussianFit);

figure(5);
heatmap(croppedImage);

%% Output Data into a .csv File
% Note that files being opened for writing should be opened with 'w'
% permission, which will delete the previous contents, or with 'a'
% permission, which will append new text to the previous contents.
fileID = fopen([pwd, '/SpotDetectorOutput.csv'],'w');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'spotID', 'plane', 'x', 'y', 'peak', 'x0', 'xdev', 'y0', 'ydev', 'ecc', 'int', 'bkgrnd');
fprintf(fileID, '%u,%u,%u,%u,%5.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.3f,%5.2f,%5.2f', spotID, planeNumber, x_coord, y_coord, ahat(1), ahat(2), ahat(3), ahat(4), ahat(5), eccentricity, unintIntensity, backgroundIntensity);
fclose(fileID);