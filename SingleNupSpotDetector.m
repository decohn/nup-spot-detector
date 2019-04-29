%% SingleNupSpotDetector
% See Github repository for README and more details regarding the project.
% Immediate term goal: continue optimizing the thresholding algorithm.
% Check the correlation between bounding box size, background intensity and total intensity,
% because I'm reasonably certain that those are very strong correlations at
% the moment. Finally, remove coordinates from my list if they are
% duplicates that point to a spot that is already pointed to.

% Go through each proposed coordinate to try to understand why some of them
% are working substantially better than others, and figure out if there's 
% anything I can do to improve the fixed bounding box condition. Then, 
% input more coordinates (both for Nup59 and for Nup60). Longer term, I
% would like to attempt a radial box and fit, but that seems pretty tough.

% On another note, it would be good here to incorporate functionality for
% the code to run on multiple proteins/image types at once. A bit of a
% hassle, but much less of one than having to switch the spot coordinates
% each time I want to switch proteins.

clc
clear variables

% excess bounding box radius
ebbr = 2;

% Images used with this program must be stored in the MATLAB home folder.
% This is easily modifiable if desired.

%% Configuration Variables 

fileID = fopen([pwd, '/SNSDconfig.txt'],'r');

imageNamePrecursor = fgetl(fileID);
spacesLocatedAt = find(imageNamePrecursor == ' ');
imageNames = imageNamePrecursor(spacesLocatedAt(2) + 1 : size(imageNamePrecursor, 2));
imageNames = split(imageNames, ",");

outputPrecursor = fgetl(fileID);
spacesLocatedAt = find(outputPrecursor == ' ');
outputFileName = outputPrecursor(spacesLocatedAt(2) + 1 : size(outputPrecursor, 2));

verbosePrecursor = fgetl(fileID);
spacesLocatedAt = find(verbosePrecursor == ' ');
verbose = string(verbosePrecursor(spacesLocatedAt(2) + 1 : size(verbosePrecursor, 2)));

% If a fixed box size is being used, then the size must be an odd perfect square
% (ex. 49, 81, etc). 9x9 seems to work well so far, because it incorporates
% a good chunk of the background.
boxSizePrecursor = fgetl(fileID);
spacesLocatedAt = find(boxSizePrecursor == ' ');
boxAreaBoundsPrecursor = boxSizePrecursor(spacesLocatedAt(2) + 1 : size(boxSizePrecursor, 2));
boxAreaBoundsPrecursor2 = split(boxAreaBoundsPrecursor, ",");
boxAreaBounds(1) = str2double(boxAreaBoundsPrecursor2{1,1});
boxAreaBounds(2) = str2double(boxAreaBoundsPrecursor2{2,1});

loopCounter = 1;
terminate = false;
while (terminate == false)
    spotCoordinatesPrecursor = fgetl(fileID);
    
    if(spotCoordinatesPrecursor ~= -1)
        spotCoordinatesPrecursor = string(spotCoordinatesPrecursor);
        spotCoordinates = split(spotCoordinatesPrecursor, ",");
        planeNumber(loopCounter) = double(spotCoordinates(1,1));
        
        % The plus one is there to adjust for the fact that arrays in
        % MATLAB don't start at the same values as FIJI arrays!
        x_coord(loopCounter) = double(spotCoordinates(2,1)) + 1;
        y_coord(loopCounter) = double(spotCoordinates(3,1)) + 1;
        imageNumber(loopCounter) = double(spotCoordinates(4,1));
        
        loopCounter = loopCounter + 1;
    else
        terminate = true;
    end
end

usableIntensities = zeros(1,loopCounter);

%% Prepare Output .csv File
% Note that files being opened for writing should be opened with 'w'
% permission, which will delete the previous contents, or with 'a'
% permission, which will append new text to the previous contents.
fileID = fopen([pwd, '/', outputFileName,'.csv'],'w');
fprintf(fileID, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'spotID', 'plane', 'x', 'y', 'peak', 'x0', 'xdev', 'y0', 'ydev', 'ecc', 'int', 'bkgrnd', 'usable', 'image number', 'bounding box area');

%% Spot Identification
for i=1:(loopCounter-1)
    
    % Image Read
    I = imread([pwd , '/', imageNames{imageNumber(i)} , '.tif'], planeNumber(i));
    
    % Scales a copy of the image in order to provide a visual reference. Should
    % only be used during testing. 
    if(verbose == "very")
        I_scaled = imadjust(I);
        figure(2*i-1);
        subplot(1,3,1);
        imshow(I_scaled, 'InitialMagnification', 'fit');
        hold on
        plot(x_coord(i), y_coord(i), 'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 20);
        title('1a');
    end

    %% Crop Image Around Identified Point
    % This code will binarize the entire image into a
    % large number of discrete objects. It will then select the single
    % object that contains the spot of interest, and display a visual reference. 
   
    if(boxAreaBounds(1) ~= boxAreaBounds(2))
        % The initial guess for the thresholding sensitivity. Setting
        % currentBoxArea in this way ensures that the loop will run at least
        % once!
        priorSensitivityGuess = 0.02;
        sensitivityGuess = 0.03;
        currentBoxArea = boxAreaBounds(2) + 1;
        boxIterations = 0;
        maxBoxIterations = 100;
        sensitivityIncrement = 0.01;

        % This loop will continuously adjust the sensitivityGuess until it
        % yields a bounding box of appropriate size (as defined by boxAreaBounds). 
        % If such a bounding box cannot be generated by any
        % value of sensitivity, the loop will give up after a given number of iterations.
        while((currentBoxArea > boxAreaBounds(2) || currentBoxArea < boxAreaBounds(1)) & boxIterations < maxBoxIterations)

            boxIterations = boxIterations + 1;
            tempSens = sensitivityGuess;

            if(currentBoxArea > boxAreaBounds(2))
                if(priorSensitivityGuess >= sensitivityGuess)
                    sensitivityGuess = sensitivityGuess - sensitivityIncrement;
                else
                    sensitivityGuess = mean([sensitivityGuess, priorSensitivityGuess]);
                end
            else
                if(priorSensitivityGuess <= sensitivityGuess)
                    sensitivityGuess = sensitivityGuess + sensitivityIncrement;
                else
                    sensitivityGuess = mean([sensitivityGuess, priorSensitivityGuess]);
                end
            end

            priorSensitivityGuess = tempSens;

            adaptivelevel = adaptthresh(I, sensitivityGuess);
            BW = imbinarize(I, adaptivelevel);
            BWOI = bwselect(BW, x_coord(i), y_coord(i), 4);

            % Gives the coordinates of the upper-left corner of the bounding box,
            % followed by the x-width and the y-width. This also includes a
            % check to ensure that the selected spot has a non-zero area.

            if(sum(sum(BWOI)) > 0)
                object_data = regionprops(BWOI, I, {'BoundingBox'});
                boundingBoxData = table2array(struct2table(object_data));
                boundingBoxData(1) = cast(boundingBoxData(1), 'uint16');
                boundingBoxData(2) = cast(boundingBoxData(2), 'uint16');

                currentBoxArea = boundingBoxData(3) * boundingBoxData(4);
            else
                currentBoxArea = 0;
            end

            % Checks to make sure that the sensitivity value isn't becoming too
            % large or too high, and forces the loop to halt if that is the
            % case.
            if(priorSensitivityGuess <= (2*sensitivityIncrement) || priorSensitivityGuess >= (1 - 2*sensitivityIncrement))
                boxIterations = maxBoxIterations;
            end
        end

        croppedImage = I((boundingBoxData(2) - ebbr) : (boundingBoxData(2) + boundingBoxData(4) - 1 + ebbr), (boundingBoxData(1) - ebbr) : (boundingBoxData(1) + boundingBoxData(3) - 1 + ebbr));
        
        % Do I have these the wrong way around?
        bboxWidth = boundingBoxData(3) + 2 * ebbr;
        bboxHeight = boundingBoxData(4) + 2 * ebbr;
    else
        % Code for a fixed bounding box size. 
        
        boxSide = sqrt(boxAreaBounds(1));
        sensitivityGuess = 0.35;
        
        adaptivelevel = adaptthresh(I, sensitivityGuess);
        BW = imbinarize(I, adaptivelevel);
        BWOI = bwselect(BW, x_coord(i), y_coord(i), 4);
        
        while(sum(sum(BWOI)) == 0)
            sensitivityGuess = sensitivityGuess + 0.01;
            adaptivelevel = adaptthresh(I, sensitivityGuess);
            BW = imbinarize(I, adaptivelevel);
            BWOI = bwselect(BW, x_coord(i), y_coord(i), 4);
        end
        
        object_data = regionprops(BWOI, I, {'Centroid'});
        centroidData = table2array(struct2table(object_data));
        centroidData(1) = floor(cast(centroidData(1), 'double'));
        centroidData(2) = floor(cast(centroidData(2), 'double')); 
        
        offset = floor(boxSide / 2.0);
        
        croppedImage = I((centroidData(2) - offset : centroidData(2) + offset), centroidData(1) - offset : centroidData(1) + offset);
        
        bboxWidth = boxSide;
        bboxHeight = boxSide;
        
        useThisSpot(i) = 1;
    end
    
    if(verbose == "very")
        subplot(1,3,2);
        imshow(BW, 'InitialMagnification', 'fit');
        hold on
    end
    
    if(verbose == "very")
        subplot(1,3,3);
        imshow(imadjust(croppedImage), 'InitialMagnification', 'fit');
        hold on
    end

    doubleCroppedImage = cast(croppedImage, 'double');

    %% Fit Cropped Image to 2D Gaussian Curve
    xdata = zeros(2, bboxWidth * bboxHeight);

    for j=1:(bboxWidth*bboxHeight)
        if(rem(j, bboxHeight) ~= 0)
            xdata(1, j) = rem(j, bboxHeight);
        else
            xdata(1, j) = bboxHeight;
        end
    end

    for j=1:bboxWidth
        xdata(2, (bboxHeight * j - (bboxHeight - 1)) : (bboxHeight * j)) = j;
    end

    ydata = reshape(doubleCroppedImage, [1, (bboxWidth  *  bboxHeight)]);

    % predicted is an anonymous fitting function that lsqcurvefit will fit the
    % data to. a will be a vector with six elements: the amplitude, the x shift, the x
    % standard deviation, the y shift, the y standard deviation, and the background value. 

    predicted = @(a, xdata) a(1) * exp(-((xdata(1, :) - a(2)).^2 / (2 * (a(3).^2))) - ((xdata(2, :) - a(4)).^2) / (2 * (a(5).^2))) + a(6);

    % a0 is the first estimate of parameters that the lsqcurvefit will use.
    % Refine this?
    a0 = [double(I(y_coord(i), x_coord(i)) - 300); 4; 4; 4; 4; 300];

    % Performs the curve fitting.
    opts = optimset('Display','off');
    [ahat, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(predicted, a0, xdata, ydata, [], [], opts);
    backgroundIntensity = ahat(6);
    
    % The final fitted Gaussian function. Can be used for the integration.
    gaussianFit = reshape(predicted(ahat, xdata), [bboxHeight, bboxWidth]);

    %% Compute Properties of the 2D Gaussian Fit
    % Use regionprops to extract the size of the portion that is above the
    % threshold?
    eccentricity = sqrt(1 - min( (ahat(3).^ 2 / ahat(5).^ 2), (ahat(5).^ 2 / ahat(3).^ 2) ) );
    
    % This more primitive method simply sums up all of the pixel values that
    % exceed the background. Then it subtracts the backgroundIntensity,
    % multiplied by the number of pixels that exceed the background. Note
    % that unintIntensity, as a uint16, is capped at 2^16 - 1. This occurs
    % I think because all of the image values are uint16 type, but this
    % could cause issues in the future.
    integratableArea = gaussianFit((1+ebbr):(bboxHeight-ebbr), (1+ebbr):(bboxWidth-ebbr));
    unintIntensity = sum(sum(integratableArea)) - backgroundIntensity * (bboxHeight - 2 * ebbr) * (bboxWidth - 2 * ebbr);

    % plots the centre of the proposed Gaussian 
    %if(verbose == "very")
     %   hold on
      %  plot(ahat(4), ahat(2), 'Marker', '.', 'LineStyle', 'none', 'MarkerSize', 20);
    %end
    
    % If a bounding box of the correct size could not be generated, don't use the spot
    if(boxAreaBounds(1) ~= boxAreaBounds(2))
        useThisSpot(i) = (boxIterations < maxBoxIterations);
    end
    
    %% Filter Spots 
    
    if(eccentricity > 0.75)
        useThisSpot(i) = 0;
    end
    
    % Sanity check.
    if(backgroundIntensity < 0)
        useThisSpot(i) = 0;
    end
    
    %% Post-Filter Cleanup
    if(verbose ~= "none" && useThisSpot(i) == 1)
        figure(2*i)
        subplot(1,2,1);
        h = heatmap(gaussianFit);
        h.Title = 'b';
        h.ColorbarVisible = 'off';

        subplot(1,2,2);
        h = heatmap(croppedImage);
        h.ColorbarVisible = 'off';
    end

    fprintf(fileID, '%u,%u,%u,%u,%5.2f,%3.2f,%3.2f,%3.2f,%3.2f,%3.3f,%5.2f,%5.2f,%u,%u,%u,\n', i, planeNumber(i), x_coord(i), y_coord(i), ahat(1), ahat(2), ahat(3), ahat(4), ahat(5), eccentricity, unintIntensity, backgroundIntensity, useThisSpot(i), imageNumber(i), bboxWidth * bboxHeight);
       
    if(useThisSpot(i) == 1)
        usableIntensities(i) = unintIntensity;
    end
end

%% Wrap things up (ex. close files) and present summary of intensities
fclose(fileID);

usableIntensities(usableIntensities == 0) = [];
numUsableSpots = size(usableIntensities, 2);

figure(2*i+1)
histogram(usableIntensities, 10, 'FaceColor', 'g');
xlabel('Integrated Intensity');
ylabel('Count of Spots');
yticks(0:7);
xlim(0:7000:7000);
% fix this so that the legend is labelled programatically
% also, the graphs should be automatically titled!
legend('Nup59');

figure(2*i+2)
scatter(ones(numUsableSpots, 1) - 0.025 + 0.05 * rand(numUsableSpots, 1), usableIntensities);
xlim([0.9, 1.1]);
xticks([]);

%figure(2*i+3)
%histfit(usableIntensities);

disp('Program completed.');