close all
materialName = 'BVO';
speciesType = 'hole';
numSpecies = 2;
inputFileName = 'unwrappedTraj.dat';
numTrajRecorded = 1.00E+02;
tFinal = 1.00E-05;
timeInterval = 1.00E-09;
bohr2ang = 0.529177249;
positionArray = dlmread('unwrappedTraj.dat') * bohr2ang;
numPathStepsPerTraj = round(tFinal / timeInterval) + 1;
positionArraySize = size(positionArray);
nSpecies = positionArraySize(2) / 3;
dataArray = zeros(numPathStepsPerTraj, numTrajRecorded * nSpecies, 3);

for trajIndex = 0:numTrajRecorded-1
    headStart = trajIndex * numPathStepsPerTraj;
    for step =0:numPathStepsPerTraj-1
        stepPosition = positionArray(headStart + step + 1, :);
        for speciesIndex = 0:nSpecies-1
            dataArray(...
                step + 1, trajIndex * nSpecies + speciesIndex + 1, :) = ...
                stepPosition(speciesIndex * 3 + 1: (speciesIndex + 1) * 3);
        end
    end
end

loops = 1000;
F(loops) = struct('cdata',[],'colormap',[]);
if numSpecies > 1
    speciesTail = 's';
else
    speciesTail = '';
end
outputFileName = strcat(materialName, '_', num2str(numSpecies), ...
                        speciesType, speciesTail, '.avi');
v = VideoWriter(outputFileName);
open(v);
index = 1;
figure('visible', 'off');
xlabelstr = sprintf('x (%c)', 197);
ylabelstr = sprintf('y (%c)', 197);
zlabelstr = sprintf('z (%c)', 197);
for step = 0:numPathStepsPerTraj-1
    if mod(step, 100) == 0 && step ~= 0
        Pext = dataArray(step + 1, :, :);
        P = reshape(Pext, 200, 3)';
        t = 0.001;
        [A , C] = MinVolEllipse(P, t);
        plot3(P(1,:),P(2,:),P(3,:),'*')
        hold on
        Ellipse_plot(A,C)
        hold off
        xlabel(xlabelstr)
        ylabel(ylabelstr)
        zlabel(zlabelstr)
        figTitle = ['Diffusion of ', num2str(numSpecies), ' ', ...
                    speciesType, speciesTail, ...
                    ' depicted over 100 traj in ms-BVO'];
        title(figTitle)
        xlim([-1000, 1000]) % bvo_2hole
        ylim([-1000, 1000]) % bvo_2hole
        zlim([-1000, 1000]) % bvo_2hole
%         xlim([-700, 700]) % bvo_2elec
%         ylim([-700, 700]) % bvo_2elec
%         zlim([-700, 700]) % bvo_2elec
%         xlim([-1400, 1400]) % hem_2elec
%         ylim([-1400, 1400]) % hem_2elec
%         zlim([-1400, 1400]) % hem_2elec
        F(index) = getframe(gcf);
        writeVideo(v,F(index));
        index = index + 1;
    end
end
close(v);