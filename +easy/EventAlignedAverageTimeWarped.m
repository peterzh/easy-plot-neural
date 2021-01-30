function grammObject = EventAlignedAverageTimeWarped(x,t,events,varargin)
%   easy.EventAlignedAverageTimeWarped(x,t,events) plots peri-event averages for the
%   continuous variable x, time-warped between the given events. Always
%   does baseline subtraction prior to the first epoch.
%
%   x should be a column vector of a continuous variable value (e.g. dF/F, 
%   pupil diameter, motion energy). 
%
%   t should be a column vector of timestamps associated with the data x. 
%
%   events should be a S x 2 cell array, where S is the number of events. 
%   First column should contain a string label for that alignment. Second column 
%   should contain a column vector of times for that event. The 2nd columns
%   should match in length between each event, and events should progress in order (one after another).
%
%   The function returns the grammObject (graphical object for the
%   figure). This function requires gramm toolbox 
%
%   easy.EventAlignedAverageTimeWarped(...,Name,Value) specifies options using one or more
%   Name,Value pair arguments.
%
%     'EpochTimePrePost' - a [1x2] row vector giving the duration of time for the
%           pre epoch (prior to 1st event) and post epoch (after last event).
%     'label'   - string which will be used to label the values (e.g
%               'z-score dF/F' or 'wheel pos');
%     'window'   - [start end] time range to plot raster and psth
%     'splitBy'   - Cell array [T x 2], one row per T split types. 
%               1st column is label, 2nd column is a label for each event
%               type. If provided, this will make a separate plot for each
%               splitting type.
%     'titleText' - Title text added to the top of the figure
%     'splitByColours' - Cell array [T x 1], where each element contains a 
%               [D x 3] matrix of RGB values where D is the number of
%               splitting conditions.
%     'ylim' - optional manual ylim applied to all plots
%     'saveFigure'   - Path to file for saving figure. PDF, PNG, EPS and JPG acceptable.

defaultPrePostEpoch = [1 1];

%Validate inputs
p = inputParser;
addRequired(p,'x',@iscolumn);
addRequired(p,'t',@iscolumn);
addRequired(p,'events',@iscell);
addParameter(p,'label','',@ischar);
addParameter(p,'splitBy',[],@iscell);
addParameter(p,'EpochTimePrePost',defaultPrePostEpoch,@(x) length(x)==2);
addParameter(p,'titleText','');
addParameter(p,'splitByColours',[],@iscell);
addParameter(p,'saveFigure',[],@ischar);
parse(p,x,t,events,varargin{:})

%Check for required toolboxes
assert(~isempty(which('gramm')),'Please install the gramm toolbox https://github.com/piermorel/gramm');

%assign variables
x = p.Results.x; %Continuous variable values
t = p.Results.t; %Continuous variable timestamps
events = p.Results.events;
splitBy = p.Results.splitBy;
EpochTimePrePost = p.Results.EpochTimePrePost;

%Check that event times given all match in length
numTrials = length(events{1,2});
assert( all( cellfun(@length,events(:,2))==numTrials ) , 'events timestamps are not all the same length');

%Calculate the warpSizes based on the ratio of time-intervals between the given events
numEpoches = size(events,1)+1;
epochTimes = [ [events{1,2}]-EpochTimePrePost(1), [events{:,2}], [events{end,2}]+EpochTimePrePost(2) ];
epochDelta = nan(1,numEpoches);
for e = 1:numEpoches
    epochDelta(e) = nanmean(epochTimes(:,1+e) - epochTimes(:,e));
end
numEpochElements = round(epochDelta*400/sum(epochDelta)); %rescale to get ~400 elements

%Extract the X value at the warped timestamps
warp_timestamps = cell(numTrials,numEpoches);
for e = 1:numEpoches
    for tr = 1:numTrials
        warp_timestamps{tr,e} = linspace(epochTimes(tr,e), epochTimes(tr,e+1), numEpochElements(e));
    end
end
warp_timestamps = cell2mat(warp_timestamps);
x_warped = interp1(t,x,warp_timestamps); %Get warped X values
x_warped = x_warped - mean(x_warped(:,1:numEpochElements(1)),2); %Baseline subtract

%Generate gramm objects for plotting
if isempty(splitBy)
    grammObject = gramm('x',1:sum(numEpochElements),'y',x_warped);
    grammObject.stat_summary('setylim', true, 'type', 'sem');
    grammObject.set_names('y',p.Results.label);
else
    clear grammObject;
    for sp = 1:size(splitBy,1)
        grammObject(sp,1) = gramm('x',1:sum(numEpochElements),'y',x_warped,'color',splitBy{sp,2});
        grammObject(sp,1).stat_summary('setylim', true, 'type', 'sem');
        grammObject(sp,1).set_names('y',p.Results.label,'color',splitBy{sp,1});
        
        if ~isempty(p.Results.splitByColours) && ~isempty(p.Results.splitByColours{sp})
            numColours = size(p.Results.splitByColours{sp},1);
            grammObject(sp,1).set_color_options('map',p.Results.splitByColours{sp},'n_color',numColours,'n_lightness',1); %default colourmap
        end
    end
end

grammObject.axe_property('TickDir','out','xtick','','xcolor','none');
grammObject.set_title(p.Results.titleText);
grammObject.set_text_options('base_size',10);
grammObject.set_layout_options('legend_width',0.2);
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
grammObject.draw();

%Label epochs
for e = 1:numEpoches-1
    arrayfun(@(ax) xline(ax, sum(numEpochElements(1:e)), 'k:',  events{e,1}), [grammObject.facet_axes_handles]);
end
arrayfun(@(ax) xline(ax, 1, 'k:',  sprintf('%s %+0.1f',events{1,1},-p.Results.EpochTimePrePost(1))), [grammObject.facet_axes_handles]);
arrayfun(@(ax) xline(ax, sum(numEpochElements), 'k:',  sprintf('%s %+0.1f',events{end,1},p.Results.EpochTimePrePost(2))), [grammObject.facet_axes_handles]);
linkaxes([grammObject.facet_axes_handles],'y');

%Save figure if desired
if ~isempty(p.Results.saveFigure)
    [fpath,fname,fext] = fileparts(p.Results.saveFigure);
    grammObject.export('file_name',fname,'export_path',fpath,'file_type',fext(2:end));
end
end