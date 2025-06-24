function [finalOut] = chat_callFreshAnalysis_inputSafety_032025(filepath, ch_height, De_np, wC, thresholds)
    % User-friendly version of callFreshAnalysis6Rec with safe input handling
%ch_height=25.15; De_np =  32.221519212327443290203718429908; wC = 10; thresholds = [15,100];filepath = '/Users/taylorthomsen/Library/CloudStorage/GoogleDrive-taylorthomsen@berkeley.edu/Shared drives/Sohn Research Lab/Project Folders/Adipocytes/adipocyte troubleshooting MZ 062025/062025_adipocytes_troubleshooting_control_80mbar_9umwc_try1/preproc';
%changed wNP in frechanalysis6rec
    %% Load data
    load(filepath, 'Rf', 'yas2det', 'sampleRate');

    if nargin < 5 || isempty(thresholds)
        thresholds = [135, 1100];
        fprintf('Default thresholds set to [%3.2f, %3.2f]\n', thresholds);
    end

    % Pad dataset
    num_points_flat = 7000;
    yas2det = padarray(yas2det, [num_points_flat, 0], 'pre');

    % Find peaks
    window = movmean(yas2det, 10000);
    [~, locs] = findpeaks(window, 'MinPeakHeight', 1000, 'MinPeakDistance', 7000); %original height 1000 and distance 9000

    finalOut = ones(length(locs), 20);

    %% Analyze each detected peak
    for i = 1:length(locs)
        fprintf('\nAnalyzing cell event %d/%d...\n', i, length(locs));

        pk = locs(i);
        startidx = pk - 4000; %original 6000
        endidx = pk + 6000;%original 8000

        if endidx > length(yas2det)
            warning('End index exceeds data length. Skipping event.');
            continue;
        end

        adjustedWindow = false;

        while true
            filteredData = movmean(yas2det(startidx:endidx), 40);
            diff1 = diff(filteredData);

            figure(1);
            %set(figure(1),'Position', [2000, 500, 2000, 400]) ;
            plot(diff1);
            title(sprintf('Difference Signal for Event %d', i));
            xlabel('Index'); ylabel('Amplitude');
            hold on;
            yline(-thresholds(1), 'r', 'LineWidth', 1.5);
            yline(-thresholds(2), 'g', 'LineWidth', 1.5);
            yline(thresholds(1), 'r', 'LineWidth', 1.5);
            yline(thresholds(2), 'g', 'LineWidth', 1.5);
            hold off;

            userInputSkip = safe_input('Do you want to skip this event? (1 for Yes, 2 for No): ', @(x) ismember(x, [1, 2]));
            if userInputSkip == 1
                break;
            end

            if ~adjustedWindow
                userInputWindow = safe_input('Do you want to adjust the window? (1 for Yes, 2 for No): ', @(x) ismember(x, [1, 2]));
                if userInputWindow == 1
                    startAdjust = safe_input('Adjust start index (relative, e.g., -1000): ', @isnumeric);
                    endAdjust = safe_input('Adjust end index (relative, e.g., 1000): ', @isnumeric);

                    startidx = max(1, startidx + startAdjust);
                    endidx = min(length(yas2det), endidx + endAdjust);
                    adjustedWindow = true;
                    continue;
                end
            end

            adjustThresholds = safe_input('Do you want to adjust thresholds? (1 for Yes, 2 for No): ', @(x) ismember(x, [1, 2]));
            if adjustThresholds == 1
                thresholds(1) = safe_input('Enter new reference threshold: ', @isnumeric);
                thresholds(2) = safe_input('Enter new squeeze threshold: ', @isnumeric);
                continue;
            end

            analyzeEvent = safe_input('Do you want to analyze this event? (1 for Yes, 2 for No): ', @(x) ismember(x, [1, 2]));
            if analyzeEvent == 2
                break;
            end

            [OUT_array] = freshAnalysis6Rec(Rf, yas2det, startidx, endidx, sampleRate, ch_height, De_np, wC, thresholds, false, false);

            if size(OUT_array, 1) ~= 1
                warning('Analysis failed for event %d. Skipping.', i);
                break;
            end

            finalOut(i, :) = OUT_array;
            save('analysis_results.mat', 'finalOut');
            fprintf('Event %d analyzed and saved.\n', i);
            break;
        end

        if userInputSkip == 1
            continue;
        end
    end

    finalOut(all(finalOut == 1, 2), :) = [];
    fprintf('\nAnalysis complete. Results saved to analysis_results.mat.\n');
end

function val = safe_input(prompt, validator)
    while true
        try
            val = input(prompt);
            if isempty(val)
                error('InputEmpty', 'No input provided.');
            end
            if nargin > 1 && ~validator(val)
                error('InvalidInput', 'Input did not pass validation.');
            end
            break;
        catch ME
            fprintf('Error: %s Please try again.\n', ME.message);
        end
    end
end
