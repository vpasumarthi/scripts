function mvee(materialName, speciesType, numSpecies, numTrajRecorded, ...
    tFinal, timeInterval, numFrames, plotPosData, ellipsoidConstruct, ...
    anisotropy, highResolution, inputFileName, bohr2ang, nDim, tol)

tic
global speciesTail timeIntervalPerFrame SEC2uS

% Construct array of positional data
positionArray = dlmread(inputFileName) * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
posDataArray = zeros(numPathStepsPerTraj, numTrajRecorded * numSpecies, ...
    nDim);
for trajIndex = 0:numTrajRecorded-1
    headStart = trajIndex * numPathStepsPerTraj;
    for step =0:numPathStepsPerTraj-1
        stepPosition = positionArray(headStart + step + 1, :);
        for speciesIndex = 0:numSpecies-1
            posDataArray(step + 1, ...
                trajIndex * numSpecies + speciesIndex + 1, :) = ...
                stepPosition(speciesIndex * nDim + 1: ...
                (speciesIndex + 1) * nDim);
        end
    end
end

% Create a video file
if numSpecies > 1
    speciesTail = 's';
else
    speciesTail = '';
end
outputVideoFileName = strcat(materialName, '_', num2str(numSpecies), ...
    speciesType, speciesTail, '.avi');
videoFile = VideoWriter(outputVideoFileName);
open(videoFile);

% Compute axes limits for frame
finalPosArray = reshape(...
    posDataArray(step + 1, :, :), numTrajRecorded * numSpecies, nDim)';
axesLimits = computeFrameLimits(finalPosArray, plotPosData, ...
    ellipsoidConstruct, nDim, numTrajRecorded, numSpecies, tol);

% Generate video frame from ellipsoid analysis
frameIndex = 1;
tempFileName = 'temp';
SEC2uS = 1.00E+06;
timeIntervalPerFrame = tFinal / numFrames * SEC2uS;
numStepsPerFrame = round((numPathStepsPerTraj - 1) / numFrames);
semiAxesLengths = zeros(numFrames, nDim);
cartesianSemiAxesLengths = zeros(numFrames, nDim);
FA = zeros(numFrames, 1);
for step = 0:numPathStepsPerTraj-1
    if mod(step, numStepsPerFrame) == 0 && step ~= 0
        stepPosData = reshape(posDataArray(step + 1, :, :), ...
            numTrajRecorded * numSpecies, nDim)';
        if ellipsoidConstruct
            sumEllipseMatrix = zeros(nDim);
            sumCenter = zeros(nDim, 1);
            sumSemiAxesLengths = zeros(1, nDim);
            sumCartesianSemiAxesLengths = zeros(1, nDim);
            for trajIndex = 0:numTrajRecorded-1
                trajStepPosData = stepPosData(...
                    :, trajIndex * numSpecies + (1:numSpecies));
                [trajEllipseMatrix , trajCenter] = MinVolEllipse(...
                    trajStepPosData, tol);
                sumEllipseMatrix = sumEllipseMatrix + trajEllipseMatrix;
                sumCenter = sumCenter + trajCenter;
                [trajSemiAxesLengths, trajCartesianSemiAxesLengths] = ...
                    axesLengths(trajEllipseMatrix);
                sumSemiAxesLengths = sumSemiAxesLengths + ...
                    trajSemiAxesLengths;
                sumCartesianSemiAxesLengths = ...
                    sumCartesianSemiAxesLengths ...
                    + trajCartesianSemiAxesLengths;
            end
            ellipseMatrix = sumEllipseMatrix / numTrajRecorded;
            center = sumCenter / numTrajRecorded;
            semiAxesLengths(frameIndex, :) = sumSemiAxesLengths / ...
                numTrajRecorded;
            cartesianSemiAxesLengths(frameIndex, :) = ...
                sumCartesianSemiAxesLengths / numTrajRecorded;
        else
            [ellipseMatrix , center] = MinVolEllipse(stepPosData, tol);
            [semiAxesLengths(frameIndex, :), ...
                cartesianSemiAxesLengths(frameIndex, :)] = ...
                axesLengths(ellipseMatrix);
        end
        timeFrame = frameIndex * timeIntervalPerFrame;
        
        % Major axes properties
        [eigVec, eigVal] = eig(ellipseMatrix);
        [~, majorIndex] = min(diag(eigVal));
        majorAxesVector = eigVec(:, majorIndex) * ...
            semiAxesLengths(frameIndex, majorIndex);
        majorAxes = [center, center] + [-majorAxesVector, majorAxesVector];
        displayMajorAxes = 1;
        
        [videoFrame, FA(frameIndex)] = generateVideoFrame(...
            numSpecies, speciesType, numTrajRecorded, materialName, ...
            ellipseMatrix, center, plotPosData, stepPosData, ...
            axesLimits, timeFrame, highResolution, tempFileName, ...
            majorAxes, displayMajorAxes);
        writeVideo(videoFile, videoFrame);
        frameIndex = frameIndex + 1;
    end
end
close(videoFile);
if highResolution
    delete([tempFileName, '.png'])
end

plotTimeEvolutionSeries(semiAxesLengths, cartesianSemiAxesLengths, ...
    anisotropy, FA, tFinal, speciesType, numTrajRecorded, numSpecies, nDim)

fileID = fopen('Run.log', 'w');
fprintf(fileID, ['Elapsed time is ', num2str(toc), ' seconds.']);
fclose(fileID);

end

function [semiAxesLengths, cartesianSemiAxesLengths] = ...
    axesLengths(ellipseMatrix)

nDim = length(ellipseMatrix);
semiAxesLengths = zeros(1, nDim);
cartesianSemiAxesLengths = zeros(1, nDim);
[eigVec, eigValMatrix] = eig(ellipseMatrix);
eigVal = diag(eigValMatrix);
semiAxesLengths(1, :) = eigVal.^-0.5;
cartesianSemiAxesLengths(1, :) = abs(eigVec * semiAxesLengths');

end

function [F, FA] = generateVideoFrame(...
    numSpecies, speciesType, numTrajRecorded, materialName, ...
    ellipseMatrix, center, plotPosData, stepPosData, boundLimits, ...
    timeFrame, highResolution, tempFileName, majorAxes, displayMajorAxes)

global speciesTail

figure('visible', 'off');
if plotPosData
    plot3(stepPosData(1,:), stepPosData(2,:), stepPosData(3,:), '*')
    hold on
end
Ellipse_plot(ellipseMatrix, center)
if displayMajorAxes
    hold on
    plot3(majorAxes(1, :), majorAxes(2, :), majorAxes(3, :), 'r')
end
hold off

xlabel(sprintf('x (%c)', 197))
ylabel(sprintf('y (%c)', 197))
zlabel(sprintf('z (%c)', 197))
figTitle = {['Time evolution of MVEE for ', num2str(numSpecies), ' ', ...
    speciesType, speciesTail, ' in ', materialName], ...
    [' averaged over ', num2str(numTrajRecorded), ' trajectories']};
title(figTitle)
xlim(boundLimits(1, :))
ylim(boundLimits(2, :))
zlim(boundLimits(3, :))

R = ellipseMatrix / trace(ellipseMatrix);
FA = sqrt((3 - (1 / trace(R^2))) / 2);
dim = [.8 .3 .3 .3];
text = {['Time = ', num2str(timeFrame), sprintf(' %cs', 956)], ...
    ['FA = ', num2str(FA)]};
annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on');

if highResolution
    print(gcf, tempFileName, '-dpng', '-r300')
    F = imread([tempFileName , '.png']);
else
    F = getframe(gcf);
end

end

function boundLimits = computeFrameLimits(finalPosArray, plotPosData, ...
    ellipsoidConstruct, nDim, numTrajRecorded, numSpecies, tol)

if plotPosData
    minPosLimit = min(min(finalPosArray, [], 2));
    maxPosLimit = max(max(finalPosArray, [], 2));
    posLimits = max(abs(minPosLimit), abs(maxPosLimit));
    numDigits = ceil(log10(posLimits));
    boundLimitValue = ceil(posLimits / 10^(numDigits - 1)) ...
        * 10^(numDigits - 1);
    boundLimits = ones(nDim, 2) .* [-1, 1] * boundLimitValue;
else
    if ellipsoidConstruct
        sumEllipseMatrix = zeros(nDim);
        sumCenter = zeros(nDim, 1);
        for trajIndex = 0:numTrajRecorded-1
            trajStepPosData = finalPosArray(...
                :, trajIndex * numSpecies + (1:numSpecies));
            [trajEllipseMatrix , trajCenter] = MinVolEllipse(...
                trajStepPosData, tol);
            sumEllipseMatrix = sumEllipseMatrix + trajEllipseMatrix;
            sumCenter = sumCenter + trajCenter;
        end
        ellipseMatrix = sumEllipseMatrix / numTrajRecorded;
        center = sumCenter / numTrajRecorded;
    else
        [ellipseMatrix , center] = MinVolEllipse(finalPosArray, tol);
    end
    [~, cartesianSemiAxesLengths] = axesLengths(ellipseMatrix);
    maxSemiAxesLength = max(cartesianSemiAxesLengths);
    posLimits = [center - maxSemiAxesLength, center + maxSemiAxesLength];
    numDigits = ceil(log10(abs(posLimits)));
    boundLimits = sign(posLimits) .* ceil(abs(posLimits) ...
        ./ 10.^(numDigits - 1)) .* 10.^(numDigits - 1);
end

end

function plotTimeEvolutionSeries(semiAxesLengths, ...
    cartesianSemiAxesLengths, anisotropy, FA, tFinal, speciesType, ...
    numTrajRecorded, numSpecies, nDim)

global speciesTail timeIntervalPerFrame SEC2uS

timeSeries = timeIntervalPerFrame:timeIntervalPerFrame:(tFinal * SEC2uS);
% Plot time evolution of ellipsoid shape
figure('visible', 'off');
plot(timeSeries, semiAxesLengths)
xlabel(sprintf('Simulation Time (%cs)', 956))
ylabel(sprintf('Magnitude of semi-axes (%c)', 197));
title('Time evolution of ellipsoid shape')
legend('a', 'b', 'c', 'Location', 'southeast')
dim = [.15 .55 .3 .3];
text = {['Num_{', speciesType, speciesTail, '} = ', ...
    num2str(numSpecies)], ['Num_{traj} = ', num2str(numTrajRecorded)]};
annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center');
figTitle = ['EllipsoidShape_', num2str(numSpecies), speciesType, ...
    speciesTail, '.png'];
saveas(gcf, figTitle)

% Plot time evolution of bounding box limits
figure('visible', 'off');
plot(timeSeries, cartesianSemiAxesLengths)
xlabel(sprintf('Simulation Time (%cs)', 956))
ylabel(sprintf('Magnitude of cartesian projected semi-axes (%c)', 197));
title('Time evolution of cartesian projected ellipsoid shape')
legend('X', 'Y', 'Z', 'Location', 'southeast')
dim = [.15 .55 .3 .3];
text = {['Num_{', speciesType, speciesTail, '} = ', ...
    num2str(numSpecies)], ['Num_{traj} = ', num2str(numTrajRecorded)]};
annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center');
figTitle = ['CartesianEllipsoidShape_', num2str(numSpecies), ...
    speciesType, speciesTail, '.png'];
saveas(gcf, figTitle)

% Plot time evolution of anisotropy
figure('visible', 'off');
dir01 = anisotropy;
dir02 = setdiff(1:nDim, anisotropy);
axesLabels = ['a', 'b', 'c'];
if length(dir01) > 1
    str01 = [axesLabels(dir01), '-plane'];
    str02 = [axesLabels(dir02), '-direction'];
    boundingLimits01 = sum(cartesianSemiAxesLengths(:, dir01).^2, 2).^0.5;
    boundingLimits02 = cartesianSemiAxesLengths(:, dir02);
else
    str01 = [axesLabels(dir01), '-direction'];
    str02 = [axesLabels(dir02), '-plane'];
    boundingLimits01 = cartesianSemiAxesLengths(:, dir01);
    boundingLimits02 = sum(cartesianSemiAxesLengths(:, dir02).^2, 2).^0.5;
end
degreeOfAnisotropy = boundingLimits01 ./ boundingLimits02;
plot(timeSeries, degreeOfAnisotropy)
xlabel(sprintf('Simulation Time (%cs)', 956))
ylabel('Degree of anisotropy')
title(['Time evolution of anisotropy in ', str01, ' vs. ', str02])
dim = [.55 .6 .3 .3];
text = {['Num_{', speciesType, speciesTail, '} = ', ...
    num2str(numSpecies)], ['Num_{traj} = ', num2str(numTrajRecorded)]};
annotation('textbox', dim, 'String', text, 'FitBoxToText', 'on', ...
    'HorizontalAlignment', 'center');
figTitle = ['DegreeOfAnistropy_', num2str(numSpecies), speciesType, ...
    speciesTail, '.png'];
saveas(gcf, figTitle)

% Plot time evolution of fractional anisotropy
figure('visible', 'off');
plot(timeSeries, FA)
xlabel(sprintf('Simulation Time (%cs)', 956))
ylabel('Fractional anisotropy')
title('Time evolution of fractional anisotropy in ellipsoid shape')
figTitle = ['FractionalAnistropy_', num2str(numSpecies), speciesType, ...
    speciesTail, '.png'];
saveas(gcf, figTitle)

end
