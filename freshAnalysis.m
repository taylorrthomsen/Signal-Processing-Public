    function [OUT_array] = ...
    freshAnalysis(Rf, yas2det, startidx, endidx, sampleRate, ch_height, De_np, wC, thresholds, plotflag, fitflag)
%% SECTION 0: device parameters
%preprocTTspet2023 is used, and variables are transformed to the one used here

% segment layout
total_segs = 12;
num_ref_segs = 1;
num_rec_segs = 10;

% check layout
if total_segs ~= (num_ref_segs + num_rec_segs + 1)
    error('check device geometry! number of segments doesn''t match!');
end

% mask geometry [um]
L = 6250; % sNPS_ver2.1 (both) - total length of the NPS channel (from start of 1st pore to end of last pore)
npL_ref = [800]; % sNPS_ver2.1 (both)
npL_rec = [290,290,290,290,290,290,290,290,290,290]; % sNPS_ver2.1 (both)
sqL = 2000; % sNPS_ver2.1 (both)
wNP = 20; % sNPS_ver2.1 (both)

% calculate De_c based on De_np
De_c = De_np*(wC/wNP)^(0.5); % D_e (effective diameter) for contraction segment

%for adipocytes, copy the folowing into command line:
%fitflag = true;
% plotflag = true;
% >> data_vector = data;
% >> ch_height = 30; De_np = 15.97;
% wC = 7;
%thresholds = [enter ref thresh, contraction thresh]
%%

%preprocess data using sept code
%look at plot and find beginning and end indeces of the cell event you want
%to analyze

%after looking at the magnitude of resistance, choose ref/rec and squeeze
%thresholds and add thresholds = [   ,   ] to workspace

Fs = sampleRate/1000; %convert to kHz
%% 
%do mov mean on cell event

% filterRf=movmean(Rf, 80); 

%choose cell event index
%Rf= Rf(startidx:endidx); 

%filteredData=movmean(yas2det, 80); 

filteredData = movmean(yas2det(startidx:endidx), 40);

%difference signal of filtered and indexed data
diff1Event = diff(filteredData); 
%% diff

% threshold values below thresholds(1)
diff1Event(abs(diff1Event) < thresholds(1)) = 0; 

% zero out the first few values
diff1Event(1:10) = 0; 

% zero out the last few values
diff1Event(end-10:end) = 0; 

% ensure all values in squeeze channel are zero
% find where squeeze starts
squeeze_begin = find(diff1Event >= thresholds(2),1); 

while (diff1Event(squeeze_begin) >= thresholds(2))
    squeeze_begin = squeeze_begin + 1;
end

% find end of squeeze
squeeze_end = find(diff1Event(squeeze_begin+20:end) <= -thresholds(2)); 

% zero out the squeeze channel
diff1Event(squeeze_begin+1:squeeze_begin+squeeze_end+20) = 0; 
%% SECTION 3: identify nonzero differences (nz_mat is the matrix of nonzero differences)

nz_mat = ones(2,length(nonzeros(diff1Event))); % preallocation
k = 1; % index for A
for i = 1:length(diff1Event)
    if (diff1Event(i) ~= 0) % look for nonzeros
        nz_mat(1,k) = i; % array index of nonzero
        nz_mat(2,k) = diff1Event(i); % nonzero value
        k = k+1;
    end
end
%% SECTION 4: remove error from A
 
k=1;
while (k < ceil(log(length(nz_mat))))
    i=1;
    while (i < length(nz_mat))

        % Case 1: current and next both positive && next > current
        if nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                nz_mat(2,i+1) > nz_mat(2,i)
            % move next into current
            nz_mat(:,i) = nz_mat(:,i+1);

        % Case 2: current and next both positive && current > next
        elseif nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                nz_mat(2,i) > nz_mat(2,i+1)
            % move current into next
            nz_mat(:,i+1) = nz_mat(:,i);

        % Case 3: current and next both negative && current > next
        elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                nz_mat(2,i) > nz_mat(2,i+1)
            % move next into current
            nz_mat(:,i) = nz_mat(:,i+1);

        % Case 4: current and next both negative && next > current
        elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                nz_mat(2,i+1) > nz_mat(2,i)
            % move current into next
            nz_mat(:,i+1) = nz_mat(:,i);
        end

        i=i+1;
    end
    k=k+1;
end
%% SECTION 5: remove repeats in A

unique_xs = unique(nz_mat(1,:)); % unique values
unique_is = ones(1,length(unique_xs)); % unique indices
for i = 1:length(unique_is)
    unique_is(i) = find(nz_mat(1,:) == unique_xs(i), 1);
end
unique_ys = nz_mat(2,unique_is);
nz_mat = [unique_xs; unique_ys];
 %% SECTION 6: rectangularize pulses

ym_rect = filteredData;
k = 1;
while (k <= 50)
    i = 1;
    while (i < length(nz_mat))

        if nz_mat(2,i) > 0 && nz_mat(2,i+1) < 0 % look for sign change in differences
            % replace all values in between with mean
            ym_rect(nz_mat(1,i):nz_mat(1,i+1)) = ...
                mean(filteredData(nz_mat(1,i):nz_mat(1,i+1))); 
        end

        i=i+1;
    end
    k=k+1;
end

  
   %% SECTION 7: Plot figures if flag is true

   % if youve already parsed out your event, max(diff1Event) will always be
   % >thresholds(1); ie: empty = 0 (false)
    if max(diff1Event) < thresholds(1) % waste of time, no pulse
        plotflag = false;
        empty = true;
    else
        empty = false;
    end

    if plotflag

        if isempty( findobj('type','figure', 'number',42) ) % create a new figure
            Pix_SS = get(0,'screensize');
            figh = figure(42);
            figsize = [0.1 0.1 0.45 0.75]*Pix_SS(4);
            set(figh,'units','pixels','pos',figsize);
        else % use the existing figure (don't change size/location)
            figh = figure(42);
        end
        
        % take top and bottom 3 values
        nsorted_d = sort(diff1Event);
        min_vals = nsorted_d(1:3);

        psorted_d = sort(diff1Event,'descend');
        max_vals = psorted_d(1:3);

        % set auto-thresholds
        if abs(min_vals(3)) < abs(max_vals(3))
            auto_thresh_value = abs(min_vals(3));
        else
            auto_thresh_value = abs(max_vals(3));
        end

        % difference plot
        ax1 = subplot(3,1,1);
        figwin_tighten();
        plot(diff1Event,'k-', 'LineWidth',1);
        title('y_{diff}');
        set(gca,'FontSize',10);
        grid(ax1,'on');
        ax1.XMinorGrid = 'on';
        axis([0, length(diff1Event), 1.1*mean(min_vals), 1.1*mean(max_vals)]);
        for i = 1:length(max_vals)
            label_str = sprintf('%3.3e',min_vals(i));
            text(i*600,1.35*max_vals(1),label_str,'FontSize',12);

            label_str = sprintf('%3.3e',max_vals(i));
            text(i*600,1.85*max_vals(1),label_str,'FontSize',12);
        end

        % plot thresholds
        hold(ax1,'on');
        for yi=[1,-1]
            yline(ax1, yi*thresholds(1), 'b', 'linew',1);
            yline(ax1, yi*thresholds(2), 'b--', 'linew',1.5);
        end
        hold(ax1,'off');

        % rectangularized
        ax2 = subplot(3,1,2);
        figwin_tighten();
        plot(ym_rect, 'k', 'linew',1);
        title('y_{rect}');
        set(gca,'FontSize',10);
        grid(ax2,'on');
        ax2.XMinorGrid = 'on';
        axis([0, length(ym_rect), 1.1*min(ym_rect), 0.01]);

        % smoothed
        ax3 = subplot(3,1,3);
        figwin_tighten();
        plot(filterRf, 'k', 'linew',1);
        hold(ax3,'on'); plot(filterRf - filteredData, 'm', 'linew',1); hold(ax3,'off'); %filterRf is "ym" and filteredData is "y_detrend"
        title('y_{LP}');
        set(gca,'FontSize',10);
        grid(ax3,'on');
        ax3.XMinorGrid = 'on';
        axis([0, length(filterRf), 0.999*min(filterRf), 1.001*max(filterRf)]);

        linkaxes([ax1,ax2,ax3], 'x');

    else
        auto_thresh_value = [];
    end

 %% SECTION 8: Detect NPS pulses
 
 %pulse_series is a matrix with the indices and parameters for rectangular pulses

    i=1;
    k = 0;
    backset = 10;
    pulse_series = ones(length(nz_mat),5);
    while (i < length(nz_mat))
        if nz_mat(2,i) > 0 && nz_mat(2,i+1) < 0 % starts positive and flips sign
            k = k + 1;
            pulse_series(k,1) = nz_mat(1,i); % Start index
            pulse_series(k,2) = nz_mat(1,i+1); % End index
            pulse_series(k,3) = mean(Rf((nz_mat(1,i)-backset):(nz_mat(1,i)-backset+10))); % normalized baseline resistance
            pulse_series(k,4) = mean(filteredData(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % avg resistance drop between pulses
            pulse_series(k,5) = std(filteredData(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % std dev of resistance drop
        end
        i=i+1;
    end

    % remove empty entries in P
    cci = 1;
    stopc = size(pulse_series,1);
    while(cci <= stopc)
        if (pulse_series(cci,1) == 1 && pulse_series(cci,2) == 1)
            pulse_series(cci,:) = [];
            stopc = stopc - 1;
        else
            cci = cci + 1;
        end
    end

    % plot detected pulses
    if plotflag
        nstarts = pulse_series(:,1);
        nstops = pulse_series(:,2);
        hold([ax1,ax2,ax3], 'all');
    
        % on the diff plot
        plot(ax1, nstarts,diff1Event(nstarts), 'bo', 'linew',1.5, 'markersi',8);
        plot(ax1, nstops,diff1Event(nstops), 'ro', 'linew',1.5, 'markersi',8);

        % on the rect plot
        for ii=1:length(nstarts)
            sc = plot(ax2, nstarts(ii):nstops(ii), ym_rect(nstarts(ii):nstops(ii)), 'linew',3);
            if mod(ii,2)==0
                set(sc, 'color','r');
            else
                set(sc, 'color','b');
            end
        end
        
        % on the data plot
        for ii=1:length(nstarts)
            sc = plot(ax3, nstarts(ii):nstops(ii), Rf(nstarts(ii):nstops(ii)), 'linew',2.5);
            if mod(ii,2)==0
                set(sc, 'color','r');
            else
                set(sc, 'color','b');
            end
        end
        
        hold([ax1,ax2,ax3], 'off');
    end

  %% SECTION 9: Extract mNPS pulse data

    % preallocate array of NPS pulse data
    out = nan(length(pulse_series) - (total_segs-1), 12);

    % preallocate array of reference segment dT_np_segs
    dT_np_segs = nan(size(out,1), num_ref_segs);

    for k = 1 : length(pulse_series)+1-total_segs

        % get relevant row numbers
        sq_k = k + num_ref_segs; % contraction segment
        ref_k_start = k; % first reference segment
        ref_k_end = sq_k - 1; % last reference segment
        rec_k_start = sq_k + 1; % first recovery segment
        rec_k_end = rec_k_start + num_rec_segs - 1; % last recovery segment
        % check segment indices
        if rec_k_end ~= (k + total_segs - 1)
            error('check segment indexing!');
        end

        %start_index = pulse_series(k,1); % starting index
        R_baseline = pulse_series(k,3); % baseline resistance
        
        % average dI & dT in reference segments

        % average node-pore Resistance drop in reference segments
        dR_np = mean(pulse_series(ref_k_start:ref_k_end,4)); %made not negative

        % average node-pore transit time in reference segments [ms]
        %   *** not applicable in JOVE device designs (sNPS_ver2.1) bc the segments are of unequal lengths
        dT_np = nan;
        % node-pore transit time in each reference segment [ms] (row vector)
        dT_nps = ( pulse_series(ref_k_start:ref_k_end,2) - pulse_series(ref_k_start:ref_k_end,1) )' ./Fs;

        % dR & dT in contraction (sqeeze) segment
        dR_c = pulse_series(sq_k,4); % squeeze resistance drop       made not negative
        dR_c_std = pulse_series(sq_k,5); % std. dev. of squeeze resistance drop
        dT_c = (pulse_series(sq_k,2) - pulse_series(sq_k,1)) /Fs; % squeeze transit time (ms)

        % post-squeeze NP resistance drops
        dR_rec = pulse_series(rec_k_start:rec_k_end, 4); % made not negative
        
      %% determine recovery time & category
        % recovery time is determined when post-squeeze NP resistance drop
        %   reaches pre-squeeze NP resistance drop (within 8% error threshold)
        % if the cell has "instant" recovery, recovery time is defined as
        %   elapsed time between the end of the squeeze segment and the
        %   beginning of the first recovery segment
        % if the cell has "transient" recovery, recovery time is defined as
        %   elapsed time between the end of the squeeze segment and the
        %   beginning of the first segment where the cell was recovered
        
        % "recovered" is when dR_rec comes within 8% of dR_np or higher
        rec_tol = 0.08;
        
        if num_rec_segs ~= 10
            warning('recovery time & category are hard-coded for devices with 10 recovery segments!');
        end
        
        rec_cat_description = '0 = instant, 1-9 = transient, 10 = prolonged';
        
        % cell was already recovered by the first recovery segment
        if (dR_np-dR_rec(1))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 0;

        % cell didn't recover until the second recovery segment
        elseif (dR_np-dR_rec(2))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+1,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 1;

        % cell didn't recover until the third recovery segment
        elseif (dR_np-dR_rec(3))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+2,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 2;

        % cell didn't recover until the fourth recovery segment
        elseif (dR_np-dR_rec(4))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+3,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 3;

        % cell didn't recover until the fifth recovery segment
        elseif (dR_np-dR_rec(5))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+4,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 4;

         % cell didn't recover until the sixth recovery segment
        elseif (dR_np-dR_rec(6))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+5,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 5;

         % cell didn't recover until the seventh recovery segment
        elseif (dR_np-dR_rec(7))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+6,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 6;

         % cell didn't recover until the eigth recovery segment
        elseif (dR_np-dR_rec(8))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+7,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 7;

          % cell didn't recover until the ninth recovery segment
        elseif (dR_np-dR_rec(9))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+8,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 8;
            
           % cell didn't recover until the tenth recovery segment
        elseif (dR_np-dR_rec(10))/dR_np < rec_tol
            T_rec = (pulse_series(rec_k_start+9,1) - pulse_series(sq_k,2)) /Fs;
            rec_cat = 9;

        % by the tenth recovery segment, cell still hadn't recovered
        else
            T_rec = Inf;
            rec_cat = 10;

        end

        %% perform mNPS-r recovery curve fitting
        if fitflag

            if num_rec_segs ~= 10
                error('mNPS-r recovery fitting is hard-coded for devices with exactly 10 recovery segments');
            end

            % populate vector rT with time-points of recovery pulses
            rT = [ 0, ...
                (pulse_series(rec_k_start,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+1,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+2,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+3,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+4,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+5,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+6,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+7,2) - pulse_series(sq_k,1)), ...
                (pulse_series(rec_k_start+8,2) - pulse_series(sq_k,1)) ] ...
                ./ Fs; % in ms

            % populate vector rdI with current drop-amplitudes of recovery
            % pulses. Anticipate that rdI should increase
            %
            % NOTE: cannot approximate size for ellipsoid particle
            rdR = dR_rec';

            % use MATLAB fit() to fit a linear polynomial to recovery data.
            % fo has fields pertaining to fit parameters:
            %   fo.p1 is the slope of the line, not sign-bounded
            %   fo.p2 is the y-intercept, not sign-bounded
            % gof has fields pertaining to goodness of fit, but only
            % gof.rsquare will be used
            [fo, gof] = fit(rT',rdR','poly1');

        else
            fo.p1 = 0;
            fo.p2 = 0;
            gof.rsquare = 0;
        end

        out(k,:) = [startidx, R_baseline, dR_np, dR_c, dR_c_std, dT_nps, dT_c, T_rec, fo.p1, fo.p2, gof.rsquare,rec_cat]; 
        out_cols = {'start_ix', 'R_baseline', 'dR_np', 'dR_c', 'dI_c_std', ...
            'dT_nps', 'dT_c', 'T_rec', 'fo_p1', 'fo_p2', 'gof_rsquare', 'rec_cat'};
        out_units = {'index', 'data units', 'data units', 'data units', 'data units', ...
            'ms',    'ms',   'ms',    '',      '',      '',        'categorical'};

        % save transit time in each reference segment
        dT_np_segs(k,:) = dT_nps;

    end

    if fitflag

        if isempty( findobj('type','figure', 'number',43) ) % create a new figure
            Pix_SS = get(0,'screensize');
            figf = figure(43);
            figsize = [0.70 0.6 1/3 2/9]*Pix_SS(4);
            set(figf,'units','pixels','pos',figsize);
        else % use the existing figure (don't change size/location)
            figf = figure(43);
        end

        title('Curve Fit'),
        ffunc = @(x) fo.p2+fo.p1*x;
        fplot(ffunc,[rT(1), rT(end)]),
        hold on, scatter(rT,rdR), hold off,

    end
    %% Section 9b: compute derived values

    calculated = zeros(size(out,1),9);

    % diameter (um)
    calculated(:,1) = ((out(:,3)./out(:,2)*De_np^2*L)./ ...
        (1+0.8*L/De_np*out(:,3)./out(:,2))).^(1/3);

    % strain (dimensionless)
    calculated(:,2) = (calculated(:,1)-wC)./calculated(:,1);

    % np velocity (mm/s = µm/ms)
    calculated(:,3) = mean(npL_ref./dT_np_segs, 2);

    % sq velocity (mm/s = µm/ms)
    calculated(:,4) = sqL./out(:,7);

    % deformed diameter (um)
    calculated(:,5) = 0.01*(pi/4*wC)*(((out(:,4)./out(:,2)*De_c^2*L)./ ...
        (1+0.8*L/De_c*out(:,4)./out(:,2))).^(1/3)).^2;

    % wCDI (dimensionless)
    calculated(:,6) = calculated(:,4)./calculated(:,3) .* calculated(:,1)./ch_height;
    

    % recovery time (ms)
    calculated(:,7) = out(:,8);

    % recovery rate
    calculated(:,8) = 0; %giving problems, set to 0, was " fo.p1 "

    % recovery category
    calculated(:,9) = out(:,12);
    


    % column names & units
    calc_cols = {'diameter', 'strain', 'V_np',          'V_c', 'def_diameter', ...
        'wCDI',         'rec_time', 'rec_rate', 'rec_cat'};
    calc_units = {'µm', 'dimensionless', 'mm/s = µm/ms', 'mm/s = µm/ms', 'µm', ...
        'dimensionless', 'ms',      '',         'categorical'};

    %% concatenate pulse data with calculated values

    OUT_array = [out(:,1:11), calculated]; % don't include rec_cat twice
    
    % set names & units for the columns
    %column_names = [out_cols(1:11), calc_cols];
    %column_units = [out_units(1:11), calc_units];
    

end
 

