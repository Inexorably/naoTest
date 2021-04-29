%% Plot historical NAO data.

% Loop through the various mutation constants historical datas
% and save the data to c0, c1, ... etc.

historicalFilenames = {'C:\Users\Shawn\Documents\GitHub\naoTest\controllers\nao_test\pops\c0\historicalData.csv' ...
    'C:\Users\Shawn\Documents\GitHub\naoTest\controllers\nao_test\pops\c1\historicalData.csv' ...
    'C:\Users\Shawn\Documents\GitHub\naoTest\controllers\nao_test\pops\c2\historicalData.csv' ...
    'C:\Users\Shawn\Documents\GitHub\naoTest\controllers\nao_test\pops\c3\historicalData.csv'};

% How many generations to look at.
N = 1000;

for ii=0:length(historicalFilenames)-1

    % Housekeeping.
    close all

    % Some options for what to plot.
    PLOT_STABLE_TIME = false;
    PLOT_ZMP_DISTANCE = true;
    PLOT_TRANS_X_DISTANCE = false;
    PLOT_COM_VELOCITY = true;
    PLOT_STDDEV_AND_MUTATION = true;
    PLOT_FITNESS_SCORE = true;

    % Load the historical data from the csv file in the default location.
    M = table2array(readtable(historicalFilenames{ii+1}));
    
    % Check out the first N entries.
    M = M(1:N,:);

    % Per code, columns are:
    %{ 
        // Column 1: Generation number
        // Column 2: Mean stable time
        // Column 3: Max stable time
        // Column 4: Mean zmp distance
        // Column 5: Max zmp distance
        // Column 6: Mean trans x distance
        // Column 7: Max trans x distance
        // Column 8: Mean COM velocity
        // Column 9: Max COM velocity
        // Column 10: Generation stddev
        // Column 11: Mutation chance
        // Column 12: Fitness score mean
        // Column 13: Fitness score max
    %}

    % Split the table into easier to read names.
    gens = M(:,1);
    timeMean = M(:,2);
    timeMax = M(:,3);
    zmpMean = M(:,4);
    zmpMin = M(:,5);
    xMean = M(:,6);
    xMax = M(:,7);
    comVMean = M(:,8);
    comVMin = M(:,9);
    stddev = M(:,10);
    mutationChance = M(:,11);
    fitnessMean = M(:,12);
    fitnessMax = M(:,13);

    % Stick the table data into a struct for easier exporting while retaining
    % field names.
    historical.gens = gens;
    historical.timeMean = timeMean;
    historical.timeMax = timeMax;
    historical.zmpMean = zmpMean;
    historical.zmpMin = zmpMin;
    historical.xMean = xMean;
    historical.xMax = xMax;
    historical.comVMean = comVMean;
    historical.comVMin = comVMin;
    historical.stddev = stddev;
    historical.mutationChance = mutationChance;
    historical.fitnessMean = fitnessMean;
    historical.fitnessMax = fitnessMax;

    % Plot the generational stable time data.
    if (PLOT_STABLE_TIME)
        figure
        semilogx(gens, timeMean, gens, timeMax);
        title('Generational Stable Time')
        xlabel('Generations')
        ylabel('Stable Time (s)')
        legend('mean', 'max') 
    end

    % Plot the related data in a subplot form.  Create a figure for the
    % subplots.
    fig = figure;

    % Plot the generational zmp distance data.
    subplot(2,2,1);
    if (PLOT_ZMP_DISTANCE)
        semilogx(gens, zmpMean, gens, zmpMin);
        grid on
        title('Generational ZMP Distance')
        xlabel('Generations')
        ylabel('ZMP Distance (m)')
        legend('mean', 'min') 
    end

    % Plot the generational COM velocity data.
    subplot(2,2,2);
    if (PLOT_COM_VELOCITY)
        semilogx(gens, comVMean, gens, comVMin);
        grid on
        title('Generational COM Velocity')
        xlabel('Generations')
        ylabel('COM Velocity (m/s)')
        legend('mean', 'min') 
    end

    % Plot the generational mutation / stddev data.
    subplot(2,2,3);
    if (PLOT_STDDEV_AND_MUTATION)
        yyaxis left
        grid on
        semilogx(gens, stddev);
        ylabel('Population standard deviation')
        yyaxis right
        semilogx(gens, mutationChance);
        title('Generational Stddev and Mutation Chance')
        xlabel('Generations')
        ylabel('Mutation chance')
        legend('stddev', 'mutation chance') 
    end

    % Plot the generational fitness data.
    subplot(2,2,4);
    if (PLOT_FITNESS_SCORE)
        semilogx(gens, fitnessMean, gens, fitnessMax);
        grid on
        title('Generational Fitness')
        xlabel('Generations')
        ylabel('Fitness Score')
        legend('mean', 'max', 'location', 'southeast') 
    end
    
    % Title the superplot.
    sgtitle("Generational data: n = 100, c = "+num2str(ii))
    
    % Save the current figure.
    saveas(fig, "generational data n = 100, c = "+num2str(ii)+".fig")
    saveas(fig, "generational data n = 100, c = "+num2str(ii)+".jpg")

    % Close the figure before saving all variables to reduce filesize.
    close all
    clear fig
    % Save historical struct into the appropriate c0, c1 etc object.
    currentFile = "c"+num2str(ii)+".mat";
    eval("c"+num2str(ii)+"=historical;");    % Slight programming evil.  This saves c with a dynamic name.
    save(currentFile);
end
%% Compare the evolutionairy trends between different mutation constants c.

% Housekeeping.
clear all

% Load the processed generational data for each mutation constant.
load("c0")
load("c1")
load("c2")
load("c3")

% Load the processed generational data for each mutation constant into an
% array so that we can iterate through it.
c = {c0 c1 c2 c3};

% Compare the characteristics for best ZMP distance, COM velocity data,
% generational stddev, and generational fitness.

% Plot the related data in a subplot form.  Create a figure for the
% subplots.
figure

% Plot the generational zmp distance data.
subplot(2,2,1);
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.zmpMin);
    hold on
end
grid on
title('Generational minimum ZMP Distance')
xlabel('Generations')
ylabel('ZMP Distance (m)')
legend('c0', 'c1', 'c2', 'c3')
hold off

% Plot the generational COM velocity data.
subplot(2,2,2);
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.comVMin);
    hold on
end
grid on
title('Generational minimum COM Velocity')
xlabel('Generations')
ylabel('COM Velocity (m/s)')
legend('c0', 'c1', 'c2', 'c3')
hold off

% Plot the generational stddev / mutation data.
subplot(2,2,3);
grid on
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.stddev);
    hold on
end
ylabel('Population standard deviation')
title('Generational Stddev')
xlabel('Generations')
legend('c0', 'c1', 'c2', 'c3')
hold off

% Plot the generational fitness data.
subplot(2,2,4);
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.fitnessMax);
    hold on
end
grid on
title('Generational max Fitness')
xlabel('Generations')
ylabel('Fitness Score')
legend('c0', 'c1', 'c2', 'c3', 'location', 'southeast')
hold off

% Title the overall figure of subplots (super title).
sgtitle('Generational data across mutation constants')

% Compare the characteristics for mean ZMP distance, COM velocity data,
% generational stddev, and generational fitness.

% Plot the related data in a subplot form.  Create a figure for the
% subplots.
figure

% Plot the generational zmp distance data.
subplot(2,2,1);
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.zmpMean);
    hold on
end
grid on
title('Generational mean ZMP Distance')
xlabel('Generations')
ylabel('ZMP Distance (m)')
legend('c0', 'c1', 'c2', 'c3')
hold off

% Plot the generational COM velocity data.
subplot(2,2,2);
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.comVMean);
    hold on
end
grid on
title('Generational mean COM Velocity')
xlabel('Generations')
ylabel('COM Velocity (m/s)')
legend('c0', 'c1', 'c2', 'c3')
hold off

% Plot the generational stddev / mutation data.
subplot(2,2,3);
grid on
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.mutationChance);
    hold on
end
title('Generational Mutation Chance')
xlabel('Generations')
ylabel('Mutation chance')
legend('c0', 'c1', 'c2', 'c3')
hold off

% Plot the generational fitness data.
subplot(2,2,4);
for (ii=1:4)
    semilogx(c{ii}.gens, c{ii}.fitnessMean);
    hold on
end
grid on
title('Generational mean Fitness')
xlabel('Generations')
ylabel('Fitness Score')
legend('c0', 'c1', 'c2', 'c3', 'location', 'southeast')
hold off

% Title the overall figure of subplots (super title).
sgtitle('Generational data across mutation constants')