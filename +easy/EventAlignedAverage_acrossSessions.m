function EventAlignedAverage_acrossSessions(mouseName,x,t,events,varargin)
%   easy.EventAlignedAverage_acrossSessions(mouseName,x,t,events) plots session-averaged 
%   peri-event averages for the continuous variable x, aligned to specific events. 
%
%   mouseName should be a [num_sessions x 1] cell array of mouse names, used
%   for averaging across sessions within each mouse, and then averaging across mice
%
%   x should be a [num_sessions x 1] cell array, where each element in the
%   array contains a continuous variable value (e.g. dF/F, pupil diameter)
%   for that session.
%
%   t should be a [num_sessions x 1] cell array containing the timestamps for x
%   events should be a [num_sessions x 1] cell array, where each element
%   contains the [num_alignments x 2] cell array outlining the event labels
%   and event times you want to align to. This [num_alignments x 2] cell
%   array is the exact same one you use to call easy.EventAlignedAverage().
%   This function requires gramm toolbox 
%
%   easy.EventAlignedAverage_acrossSessions(...,Name,Value) specifies options using one or more
%   Name,Value pair arguments.
%
%     'label'   - string which will be used to label the values (e.g
%               'z-score dF/F' or 'wheel pos');
%     'window'   - [start end] time range to plot raster and psth
%     'baselineSubtract' - [start end] time range of baseline subtraction
%                           window
%     'splitBy'   - [num_sessions x 1] cell array, and each element is a 
%               [S x 2] cell array defining the splitting conditions. Same as 
%               the splitBy definition used in easy.EventAlignedAverage();
%     'titleText' - Title text added to the top of the figure
%     'splitByColours' - [S x 1] cell array, where each contains a 
%               [D x 3] matrix of RGB values where D is the number of
%               splitting conditions across all sessions.
%     'ylim' - optional manual ylim applied to all figures
%     'saveFigure'   - Path to file for saving figure. PDF, PNG, EPS and JPG acceptable.

default_window = [-0.5 +0.5];

%Validate inputs
p = inputParser;
addRequired(p,'mouseName',@iscell);
addRequired(p,'x',@iscell);
addRequired(p,'t',@iscell);
addRequired(p,'events',@iscell);
addParameter(p,'label','',@ischar);
addParameter(p,'splitBy',[],@iscell);
addParameter(p,'window',default_window,@isrow);
addParameter(p,'baselineSubtract',[],@(x) size(x,1)==1 & size(x,2)==2 & diff(x)>0);
addParameter(p,'titleText','');
addParameter(p,'ylim',[],@(x) length(x)==2 & diff(x)>0);
addParameter(p,'splitByColours',[],@iscell);
addParameter(p,'saveFigure',[],@ischar);
parse(p,mouseName,x,t,events,varargin{:})

%Check for required toolboxes
assert(~isempty(which('gramm')),'Please install the gramm toolbox https://github.com/piermorel/gramm');

%assign variables
x = p.Results.x; %Continuous variable values
t = p.Results.t; %Continuous variable timestamps
assert( all(cellfun(@iscolumn,x)) & all(cellfun(@iscolumn,t)), 'x and t should be a cell array of column vectors'); %Check for col vector
events = p.Results.events;
splitBy = p.Results.splitBy;
baselineSubtract = p.Results.baselineSubtract;
num_alignments = size(events{1},1);
num_sessions = length(x);
[mice,~,mouseID] = unique(p.Results.mouseName);
num_mice = length(mice);

%For each session and each event type, get aligned values of X
bins = p.Results.window(1):median( diff(cat(1,t{:})) ):p.Results.window(2);
aligned_x_perSession = cell(num_sessions,num_alignments);
splitGroups_perSession = cell(num_sessions,num_alignments);
mouseID_perSession = cell(num_sessions,num_alignments);
for i = 1:num_alignments
    for sess = 1:num_sessions
        %get aligned values around each event time
        aligned_x = interp1(t{sess},x{sess}, events{sess}{i,2} + bins);
        
        %baseline subtract (if required)
        if ~isempty(baselineSubtract)
            aligned_x = aligned_x - mean(aligned_x(:,baselineSubtract(1) <= bins & bins <= baselineSubtract(2)),2);
        end
        
        %Now average across trials, respecting splitBy groupings
        if isempty(splitBy) || isempty(splitBy{sess}{i,2}) %if no splitting
            aligned_x = nanmean(aligned_x,1);
            splitGroups_perSession{sess,i} = 1;
        else %if splitting by groups
            [splitGroups_perSession{sess,i},~,groupID] = unique(splitBy{sess}{i,2});
            aligned_x = grpstats(aligned_x, groupID, 'nanmean');
        end
        
        %Concatenate average to larger matrix for later plotting
        aligned_x_perSession{sess,i} = aligned_x;
        
        %add mouse ID, (useful later for plotting)
        mouseID_perSession{sess,i} = repmat(mouseID(sess), length(splitGroups_perSession{sess,i}), 1);
    end
end

%Check for any splitting conditions which don't have repeats
for i = 1:num_alignments
    counts = tabulate(string( cat(1,splitGroups_perSession{:,i})) );
    noRepeats = find([counts{:,2}]==1);
    if ~isempty(noRepeats)
        for c = 1:length(noRepeats)
            warning('For alignment %s, condition "%s=%s" has no repeats. Will not plot',events{1}{i,1}, splitBy{1}{i,1}, counts{noRepeats(c),1} );
        end
    end
end

%Now create averages across sessions within each mouse
aligned_x_perMouse = cell(num_mice,num_alignments);
splitGroups_perMouse = cell(num_mice,num_alignments);
splitGroupColours_perMouse = cell(num_mice,num_alignments);
splitGroupColours_acrossMice = cell(num_alignments);
for i = 1:num_alignments
    
    for m = 1:num_mice
        x = cat(1,aligned_x_perSession{mouseID==m,i});
        s = cat(1,splitGroups_perSession{mouseID==m,i});
        aligned_x_perMouse{m,i} = grpstats(x, string(s), 'mean');
        splitGroups_perMouse{m,i} = unique(s);
    end
    
    %get correct colour set for this splitting
    if ~isempty(p.Results.splitByColours) && ~isempty(p.Results.splitByColours{i})
        allSplits = cellfun(@(s) s{i,2}, splitBy, 'UniformOutput', false);
        allSplits = unique(cat(1,allSplits{:}));
    
        %Per mouse, match up the splitting conditions to the set of colours
        for m = 1:num_mice
            idx = ismember(allSplits,splitGroups_perMouse{m,i}); %get colours for the present conditions            
            splitGroupColours_perMouse{m,i} = p.Results.splitByColours{i}(idx,:);
        end
    end
end

%Plots: 

%average across sessions within each mouse
clear grammObject;
for subj = 1:num_mice
    for i = 1:num_alignments
        %create gramm psth object
        grammObject(subj,i) = gramm('x',bins,'y',cat(1,aligned_x_perSession{mouseID==subj,i}),'color',cat(1,splitGroups_perSession{mouseID==subj,i}));
        grammObject(subj,i).stat_summary('setylim', true, 'type', 'sem');
        
        %title for this alignment
        if subj == 1
            title = [events{sess}(i,1),' '];
            if ~isempty(splitBy) && ~isempty(splitBy{sess}{i,2})
                title{2} = ['split by ' splitBy{sess}{i,1}];
            end
            grammObject(subj,i).set_title(title);
        end
        
        if i == 1
            grammObject(subj,i).set_names('y',p.Results.label,'x','','color','');
        else
            grammObject(subj,i).set_names('y','','x','','color','');
        end
        grammObject(subj,i).axe_property('XLim',p.Results.window);
        
        %Set colour options
        if ~isempty(splitGroupColours_perMouse{subj,i})
            numColours = size(splitGroupColours_perMouse{subj,i},1);
            grammObject(subj,i).set_color_options('map',splitGroupColours_perMouse{subj,i},'n_color',numColours,'n_lightness',1); %default colourmap
        end
    end
end

%average across mice
if num_mice > 1
    for i = 1:num_alignments
        grammObject(num_mice+1,i) = gramm('x',bins,'y',cat(1,aligned_x_perMouse{:,i}),'color',cat(1,splitGroups_perMouse{:,i}) );
        grammObject(num_mice+1,i).stat_summary('setylim', true, 'type', 'sem');
        
        if i == 1
            grammObject(num_mice+1,i).set_names('y',p.Results.label,'x',events{1}{i,1},'color','');
        else
        	grammObject(num_mice+1,i).set_names('y','','x',events{1}{i,1},'color','');
        end
        
        %Set colour options
        if ~isempty(p.Results.splitByColours) && ~isempty(p.Results.splitByColours{i})
            numColours = size(p.Results.splitByColours{i},1);
            grammObject(num_mice+1,i).set_color_options('map',p.Results.splitByColours{i},'n_color',numColours,'n_lightness',1); %default colourmap
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
grammObject.draw();

%match all axes limits
linkaxes([grammObject.facet_axes_handles],'xy');

%If ylim provided, set all YLims
if ~isempty(p.Results.ylim)
    set([grammObject.facet_axes_handles],'ylim',p.Results.ylim);
end

%add mouse label to each plot
for m = 1:num_mice
    axPos = get(grammObject(m,1).facet_axes_handles,'position');
    annotation('textbox', axPos, 'string',mice{m},'HorizontalAlignment','left','EdgeColor','none','fontsize',15);
end
if num_mice > 1
axPos = get(grammObject(num_mice+1,1).facet_axes_handles,'position');
    annotation('textbox', axPos, 'string','Grand avg.','HorizontalAlignment','left','EdgeColor','none','fontsize',15);
end

%add vertical line to every axes at 0
arrayfun(@(ax) xline(ax,0,'k--'), [grammObject.facet_axes_handles]);

%Save figure if desired
if ~isempty(p.Results.saveFigure)
    [fpath,fname,fext] = fileparts(p.Results.saveFigure);
    grammObject.export('file_name',fname,'export_path',fpath,'file_type',fext(2:end));
end

end
