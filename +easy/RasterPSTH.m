function grammObject = RasterPSTH(spikeTimes, events, varargin)
%   easy.RasterPSTH(spikeTimes,events) plots a raster and psth for the
%   spike times aligned to an event. spikeTimes should be a column vector of
%   S x 1 spike time events. events should be a S x 2 cell array, where S 
%   is the number of alignments required. First column should contain a string label
%   for that alignment. Second column should contain a column vector of
%   times for that event.
%   The function returns the grammObject (graphical object for the
%   figure). This function requires gramm and spikes toolboxes to be
%   installed.
%
%   easy.RasterPSTH(...,Name,Value) specifies options using one or more
%   Name,Value pair arguments.
%
%     'window'   - [start end] time range to plot raster and psth. Default
%                   is [-1 1]
%     'splitBy'   - Cell array [S x 2], one row per split type. 
%               1st column is label, 2nd column is a label for each event type.
%     'sortBy'  - Cell array [S x 2]. 1st column is sorting label, 2nd
%               column is sorting times.
%     'kilosortTemplate' - optional matrix [numChannels x numTimebins]
%       giving the kilosort template. If provided, this will be plotted on
%       the left side.
%     'titleText' - Title text
%     'psthBinWidth' - Time bin width in seconds, used for computing the
%     PSTH
%     'psthSmoothWidth' - Time width of a causal smoothing kernel
%     'splitByColours' - Cell array [S x 1], where each element contains a 
%               [D x 3] matrix of RGB values where D is the number of
%               splitting conditions.
%     'ylim' - [min max] optional manual y limit for all axes.
%     'saveFigure'   - Path to file for saving figure. PDF, PNG, EPS and JPG acceptable.

default_window = [-1 +1];
default_psthBinWidth = 1/1000;
default_psthSmoothWidth = 50/1000;

%Validate inputs
p = inputParser;
addRequired(p,'spikeTimes',@(x) iscolumn(x) & ~isempty(x));
addRequired(p,'events',@iscell);
addParameter(p,'splitBy',[],@iscell);
addParameter(p,'sortBy',[],@iscell);
addParameter(p,'psthBinWidth',default_psthBinWidth,@isscalar);
addParameter(p,'psthSmoothWidth',default_psthSmoothWidth,@isscalar);
addParameter(p,'kilosortTemplate',[],@ismatrix);
addParameter(p,'window',default_window,@isrow);
addParameter(p,'ylim',[],@(x) length(x)==2 & diff(x)>0);
addParameter(p,'titleText','');
addParameter(p,'splitByColours',[],@iscell);
addParameter(p,'saveFigure',[],@ischar);
parse(p,spikeTimes,events,varargin{:})

%Check for required toolboxes
assert(~isempty(which('gramm')),'Please install the gramm toolbox https://github.com/piermorel/gramm');
assert(~isempty(which('psthAndBA')),'Please install the spikes toolbox https://github.com/cortex-lab/spikes/');

%Assign variables
kilosortTemplate = p.Results.kilosortTemplate;
spikeTimes = p.Results.spikeTimes;
events = p.Results.events;
num_alignments = size(p.Results.events,1);

if isempty(p.Results.splitBy)
    splitBy = cell(num_alignments,2);
else
    splitBy = p.Results.splitBy;
end

if isempty(p.Results.sortBy)
    sortBy = cell(num_alignments,2);
else
    sortBy = p.Results.sortBy;
end

%For each event type, compute the aligned spike times and binned psths
bins                    = cell(num_alignments,1);
binnedArray_smoothed    = cell(num_alignments,1);
spikeTimes_byEvent      = cell(num_alignments,1);
for i = 1:num_alignments
    [bins{i}, binnedArray_smoothed{i}, spikeTimes_byEvent{i}] = alignEventTimes(spikeTimes, events{i,2}, p.Results.window, p.Results.psthBinWidth, p.Results.psthSmoothWidth);
    
    %Trim binned array on the left side because smoothing window
    numBinsToTrim = round((p.Results.psthSmoothWidth*2)/p.Results.psthBinWidth);
    binnedArray_smoothed{i}(:,1:numBinsToTrim) = NaN;
    
    %If splitby and sortBy are defined, then rearrange binned array & spike
    %times by the sorting
    if ~isempty(sortBy{i,2})
        sortBy{i,2} = sortBy{i,2} - events{i,2}; %Calculate sort time relative to event time
        [sortBy{i,2},sortIdx] = sort(sortBy{i,2},'ascend');
        
        spikeTimes_byEvent{i} = spikeTimes_byEvent{i}(sortIdx);
        binnedArray_smoothed{i} = binnedArray_smoothed{i}(sortIdx,:);
                
        if ~isempty(splitBy{i,2})
            splitBy{i,2} = splitBy{i,2}(sortIdx);
        end
    end
end

clear grammObject;

%add placeholder for template plot
grammObject(1,1) = gramm('x',1000,'y',1000);
grammObject(1,1).set_title('Kilosort template');

%get spike autocorrelogram data
binSize = 0.5/1000; 
b = 0.0001:binSize:0.05; %bins for ACG
[ACG_counts,~] = histdiff(spikeTimes, spikeTimes, b);
ACG_counts = ACG_counts./binSize; 
ACG_counts = ACG_counts(:);
b = b(:);
ACG_xx = [b(1); reshape([b(2:end-1) b(2:end-1)]', (numel(b)-2)*2,1); b(end)];
ACG_yy = reshape([ACG_counts ACG_counts]', numel(ACG_counts)*2,1);
    
%Plot autocorrelation
grammObject(2,1) = gramm('x',1000*ACG_xx,'y',ACG_yy);
grammObject(2,1).geom_line();
grammObject(2,1).geom_vline( 'xintercept', 2 ); %add refractory time as 2 ms
grammObject(2,1).set_names('x','Time (ms)','y','');
grammObject(2,1).set_title('Autocorrelation');
grammObject(2,1).axe_property('Ytick','','Ycolor','none');

%For each event type, create gramm visual object
for i = 1:num_alignments
    
    %title for this alignment
    title = events(i,1);
    if ~isempty(splitBy{i,2})
        title = [title,['split by ' splitBy{i,1}]];
    end
     if ~isempty(sortBy{i,2})
        title = [title,['sorted by ' sortBy{i,1}]];
    end
    
    %create gramm raster object
    grammObject(1,1+i) = gramm('x',spikeTimes_byEvent{i},'color',splitBy{i,2});
    grammObject(1,1+i).geom_raster('geom','point');
    grammObject(1,1+i).geom_vline( 'xintercept', 0 );
    grammObject(1,1+i).set_names('y','Event number','x',events{i,1},'color','');
    grammObject(1,1+i).set_title(title);
    grammObject(1,1+i).set_point_options('base_size',2);
    grammObject(1,1+i).axe_property('XLim',p.Results.window,'Ydir','reverse');
    grammObject(1,1+i).axe_property('Ytick','','Ycolor','none');

    %create gramm psth object
    grammObject(2,1+i) = gramm('x',bins{i},'y',binnedArray_smoothed{i},'color',splitBy{i,2});
    grammObject(2,1+i).stat_summary('setylim', true,'type', 'sem');
    grammObject(2,1+i).geom_vline( 'xintercept', 0 );
    grammObject(2,1+i).set_names('y','Spikes/sec','x',events{i,1},'color','');
    grammObject(2,1+i).axe_property('XLim',p.Results.window);
    
        
    %Set colour options
    if ~isempty(p.Results.splitByColours) && ~isempty(p.Results.splitByColours{i})
        numSplits = length(unique(splitBy{i,2}));
        numColours = size(p.Results.splitByColours{i},1);
        assert(numSplits == numColours, sprintf('%s: %d colours specified but %d splitting conditions',splitBy{i,1},numColours,numSplits));
        grammObject(:,1+i).set_color_options('map',p.Results.splitByColours{i},'n_color',numColours,'n_lightness',1); %default colourmap
    end
end

%Render figure;
grammObject.set_title(p.Results.titleText);
grammObject.axe_property('TickDir','out');
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
grammObject.draw();

% Add sorting icons. Note in gramm this can only be done after figure
% rendering.
for i = 1:num_alignments
    if ~isempty(sortBy{i,2})
        grammObject(1,1+i).update('x',num2cell(sortBy{i,2}));
        grammObject(1,1+i).geom_raster('geom','point');
        grammObject(1,1+i).no_legend();
    end
end

if any(~cellfun(@isempty,sortBy(:,2)))
    grammObject.draw();
end

%Change icons to black triangles
grammObjectRasters = grammObject(1,2:end);
hasSort = ~cellfun(@isempty,sortBy(:,2));
if any(hasSort)
    hx = [grammObjectRasters(1,hasSort).results];
    for i = 1:length(hx)
        set([hx(i).geom_raster_handle],'Marker','v','MarkerSize',4,'MarkerFaceColor',[0 0 0]);
    end
end

%Match all PSTH firing rate axes
if num_alignments > 1
    linkaxes([grammObject(2,2:end).facet_axes_handles],'y');
end

%If ylim provided, set all YLims
if ~isempty(p.Results.ylim)
    set([grammObject(2,2:end).facet_axes_handles],'ylim',p.Results.ylim);
end

%Add template imagesc plot if provided
if ~isempty(kilosortTemplate)
    template_time = 1000*linspace(-40,41,82)/30000; %assumes 30,000 Hz sampling
    thisAx = grammObject(1,1).facet_axes_handles;
    imagesc( thisAx, template_time,[],kilosortTemplate )
    set(thisAx, 'YDir', 'normal', 'xlim', [-0.8 0.8],  'ylim', [0  size(kilosortTemplate,1)],'xcolor','none');
    ylabel(thisAx,'Channel Number (0 = tip)');
    colormap(thisAx,colormap_BlueWhiteRed); caxis(thisAx,[-1 1]*max(abs(caxis()))/2);
    
    %add small inset plot showing curve of waveform
    [~,top3channels] = sort(mean(abs(kilosortTemplate),2),'descend'); top3channels = top3channels(1:3);
    insetAx = axes('position',thisAx.Position);
    insetAx.Position(3:4) = insetAx.Position(3:4)/3;
    insetAx.Position(1:2) = insetAx.Position(1:2) + 0.015;
    plot(insetAx, template_time, mean(kilosortTemplate(top3channels,:),1) ,'k-','linewidth',2);
    set(insetAx,'Color','none','box','off','ycolor','none','ylim',[-1 1]*max(max(abs(kilosortTemplate))));
    arrayfun(@(x) xline(insetAx,x,'k:'),[-0.5 0 0.5]);
    arrayfun(@(x) xline(thisAx,x,'k:'),[-0.5 0 0.5]);
    xlabel(insetAx,'Time (ms)');
else
    delete(grammObject(1,1).facet_axes_handles);
    delete(grammObject(1).title_axe_handle);
end

%add small text indicating psth parameters
text = sprintf('%d ms binning\n%d ms smoothing',1000*p.Results.psthBinWidth, 1000*p.Results.psthSmoothWidth);
axPos = grammObject(2,2).facet_axes_handles.Position; %first psth plot handle position
annotation('textbox', axPos, 'string',text,'HorizontalAlignment','left','EdgeColor','none');

%Save figure if desired
if ~isempty(p.Results.saveFigure)
    [fpath,fname,fext] = fileparts(p.Results.saveFigure);
    grammObject.export('file_name',fname,'export_path',fpath,'file_type',fext(2:end));
end
end

function [bins, binnedArray_smoothed, spikeTimes_byEvent] = alignEventTimes(spikeTimes, eventTimes, window, psthBinWidth, psthSmoothWidth)
%This function computes two things: 1) spike times in a window surrounding
%an event time. 2) binned spike counts around the window (useful for PSTH
%later).
% Inputs:
% spikeTimes - column vector of spike times
% eventTimes - column vector of event times
% window - 1 x 2 vector of windows
% psthBinWidth - psth binning in seconds
% psthSmoothWidth - psth smoothing window in seconds

% Calculate binned array of spike counts
[~, bins, ~, ~, ~, binnedArray] = psthAndBA(spikeTimes, eventTimes, window, psthBinWidth);

%Set NaN event times to have NaN spike counts
binnedArray(isnan(eventTimes),:) = NaN;

% Smooth binned array using a causal Gaussian window
smoothFilt = myGaussWin(psthSmoothWidth, 1/psthBinWidth);
smoothFilt(1:round(numel(smoothFilt)/2)-1) = 0;  %Truncate before 0 sec, to generate causal filter
smoothFilt = smoothFilt./sum(smoothFilt);
binnedArray_smoothed = conv2(smoothFilt,1,binnedArray', 'same')';

%convert smoothed spike counts to firing rate
binnedArray_smoothed = binnedArray_smoothed./psthBinWidth;

% Convert vector of spike times into a cell array of spike times grouped by
% each event instance
out = WithinRanges(spikeTimes, eventTimes + window,(1:length(eventTimes))','matrix');
spikeTimes_byEvent = arrayfun( @(n) spikeTimes(logical(out(:,n))) - eventTimes(n) , 1:length(eventTimes), 'uni', 0)';
end