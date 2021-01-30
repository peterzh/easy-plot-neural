function grammObject = EventAlignedAverage(x,t,events,varargin)
%   easy.EventAlignedAverage(x,t,events) plots peri-event averages for the
%   continuous variable x, aligned to specific events. 
%
%   x should be a column vector of a continuous variable value (e.g. dF/F, 
%   pupil diameter, motion energy). 
%
%   t should be a column vector of timestamps associated with the data x. 
%
%   events should be a S x 2 cell array, where S is the number of alignments required. 
%   First column should contain a string label for that alignment. Second column 
%   should contain a column vector of times for that event.
%   The function returns the grammObject (graphical object for the
%   figure). This function requires gramm toolbox 
%
%   easy.EventAlignedAverage(...,Name,Value) specifies options using one or more
%   Name,Value pair arguments.
%
%     'label'   - string which will be used to label the values (e.g
%               'z-score dF/F' or 'wheel pos');
%     'window'   - [start end] time range to plot raster and psth
%     'baselineSubtract' - [start end] time range of baseline subtraction
%                           window
%     'splitBy'   - Cell array [S x 2], one row per split type. 
%               1st column is label, 2nd column is a label for each event
%               type. 
%     'titleText' - Title text added to the top of the figure
%     'mode' - Plot type: 'ci' (mean & 95% CI) , 'sem' (mean & standard
%               error), 'std' (mean and stdev), '95percentile' (median & 95%
%               percentiles)
%     'splitByColours' - Cell array [S x 1], where each element contains a 
%               [D x 3] matrix of RGB values where D is the number of
%               splitting conditions.
%     'ylim' - optional manual ylim applied to all plots
%     'saveFigure'   - Path to file for saving figure. PDF, PNG, EPS and JPG acceptable.

default_window = [-0.5 +0.5];

%Validate inputs
p = inputParser;
addRequired(p,'x',@iscolumn);
addRequired(p,'t',@iscolumn);
addRequired(p,'events',@iscell);
addParameter(p,'label','',@ischar);
addParameter(p,'splitBy',[],@iscell);
addParameter(p,'window',default_window,@isrow);
addParameter(p,'baselineSubtract',[],@(x) size(x,1)==1 & size(x,2)==2 & diff(x)>0);
addParameter(p,'titleText','');
addParameter(p,'mode','ci',@(x) ismember(x,{'ci','sem','std','quartile','95percentile'}));
addParameter(p,'ylim',[],@(x) length(x)==2 & diff(x)>0);
addParameter(p,'splitByColours',[],@iscell);
addParameter(p,'saveFigure',[],@ischar);
parse(p,x,t,events,varargin{:})

%Check for required toolboxes
assert(~isempty(which('gramm')),'Please install the gramm toolbox https://github.com/piermorel/gramm');

%assign variables
x = p.Results.x; %Continuous variable values
t = p.Results.t; %Continuous variable timestamps
events = p.Results.events;
num_alignments = size(p.Results.events,1);

if isempty(p.Results.splitBy)
    splitBy = cell(num_alignments,2);
else
    splitBy = p.Results.splitBy;
end

%For each event type, compute the aligned averages
bins = p.Results.window(1):mean(diff(t)):p.Results.window(2);
aligned_x = cell(num_alignments,1);
for i = 1:num_alignments
    
    %get x value around each event time
    aligned_x{i} = interp1(t,x,events{i,2} + bins);
    
    %if specified, perform baseline subtraction for each event
    if ~isempty(p.Results.baselineSubtract)
        baseline_idx = p.Results.baselineSubtract(1) <= bins & bins <= p.Results.baselineSubtract(2);
        
        aligned_x{i} = aligned_x{i} - mean(aligned_x{i}(:,baseline_idx),2);
    end
end

clear grammObject;

%For each event type, create gramm visual object
for i = 1:num_alignments
    
    %title for this alignment
    title = [events(i,1),' '];
    if ~isempty(splitBy{i,2})
        title{2} = ['split by ' splitBy{i,1}];
    end

    %create gramm psth object
    grammObject(1,i) = gramm('x',bins,'y',aligned_x{i},'color',splitBy{i,2});
    grammObject(1,i).stat_summary('setylim', true, 'type', p.Results.mode);
    grammObject(1,i).set_names('y',[p.Results.label sprintf(' [%s]',p.Results.mode)],'x',events{i,1},'color','');
    grammObject(1,i).set_title(title);
    grammObject(1,i).axe_property('XLim',p.Results.window);
    
    %Set colour options
    if ~isempty(p.Results.splitByColours) && ~isempty(p.Results.splitByColours{i})
        numSplits = length(unique(splitBy{i,2}));
        numColours = size(p.Results.splitByColours{i},1);
        assert(numSplits == numColours, sprintf('%s: %d colours specified but %d splitting conditions',splitBy{i,1},numColours,numSplits));
        grammObject(1,i).set_color_options('map',p.Results.splitByColours{i},'n_color',numColours,'n_lightness',1); %default colourmap
    end

end

%Render figure;
grammObject.set_title(p.Results.titleText);
grammObject.axe_property('TickDir','out');
grammObject.set_text_options('base_size',20);
figure('units','normalized','outerposition',[0 0 1 1],'color','w')
grammObject.draw();

%Match all y axes
if num_alignments > 1
    linkaxes([grammObject.facet_axes_handles],'y');
end

%If ylim provided, set all YLims
if ~isempty(p.Results.ylim)
    set([grammObject.facet_axes_handles],'ylim',p.Results.ylim);
end

%add vertical line to every axes at 0
arrayfun(@(ax) xline(ax,0,'k--'), [grammObject.facet_axes_handles]);

%Save figure if desired
if ~isempty(p.Results.saveFigure)
    [fpath,fname,fext] = fileparts(p.Results.saveFigure);
    grammObject.export('file_name',fname,'export_path',fpath,'file_type',fext(2:end));
end
end