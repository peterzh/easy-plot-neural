function PSTHbyDepth(ephysStruct, events, varargin)
%   easy.PSTHbyDepth(ephysStruct,events) plots psths by probe depth, for
%   different time alignments. ephysStruct should be the struct outputted
%   from loadKSDir function. events should be a S x 2 cell array, where S 
%   is the number of alignments required. First column should contain a string label
%   for that alignment. Second column should contain a column vector of
%   times for that event. This function requires spikes toolboxes to be
%   installed.
%
%   easy.PSTHbyDepth(...,Name,Value) specifies options using one or more
%   Name,Value pair arguments.
%
%     'window'   - [start end] time range to plot raster and psth
%     'depthBinSize' - Depth bin in units of channel coords (um usually)
%     'titleText' - Title text
%     'psthBinWidth' - Time bin width in seconds
%     'baselineNormalisatoin' - Time [start end] for normalising PSTHs by baseline
%     'saveFigure'   - Path to file for saving figure. PDF, PNG, EPS and JPG acceptable.

default_window = [-0.5 +0.5];
default_psthBinWidth = 10/1000;
default_depthBin = 80;
default_baselineNorm = [-0.2 -0.05];

%Validate inputs
p = inputParser;
addRequired(p,'ephysStruct',@isstruct);
addRequired(p,'events',@iscell);
addParameter(p,'depthBinSize',default_depthBin,@isscalar);
addParameter(p,'psthBinWidth',default_psthBinWidth,@isscalar);
addParameter(p,'window',default_window,@isrow);
addParameter(p,'baselineNorm',default_baselineNorm,@isrow);
addParameter(p,'titleText','');
addParameter(p,'saveFigure',[],@ischar);
parse(p,ephysStruct,events,varargin{:})

%Check for required toolboxes
assert(~isempty(which('psthAndBA')),'Please install the spikes toolbox https://github.com/cortex-lab/spikes/');

%Compute spike depths
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(ephysStruct.temps, ephysStruct.winv, ephysStruct.ycoords, ephysStruct.spikeTemplates, ephysStruct.tempScalingAmps);

%For each alignment
num_alignments = size(p.Results.events,1);
f = figure('units','normalized','outerposition',[0 0 1 1],'color','w');
ha = tight_subplot(1,num_alignments,0.02,[0.1 0.05],0.05);
for i = 1:num_alignments
    %Compute the normalised psths for each depth bin
    [timeBins, depthBins, allP, normVals] = psthByDepth(ephysStruct.st, spikeDepths, ...
        p.Results.depthBinSize, p.Results.psthBinWidth,...
        events{i,2}, p.Results.window, p.Results.baselineNorm);

    %plot 
    hold(ha(i),'on');
    imagesc(ha(i),timeBins, depthBins(1:end-1), allP);
    xline(ha(i),0,'k--','Linewidth',2);
    title(ha(i),events{i,1});
    
    if i == num_alignments
        h = colorbar;
        h.Label.String = 'Firing rate z-score';
    end
end

set(ha,'ydir','normal','clim',[-10 10],'xlim',p.Results.window,'box','off','ylim',[0 3800]);
ylabel(ha(1),'Position on electrode array (µm; 0 = tip)');
set(ha(2:end),'ytick','','ycolor','none');
linkaxes(ha,'xy');
xlabel(ha(1),'Time (sec)');
colormap(f, colormap_BlueWhiteRed);

%Add titles
annotation('textbox', [0, 1, 0.4,0], 'string', sprintf('%d ms time bin\n%d depth bin',1000*p.Results.psthBinWidth,p.Results.depthBinSize),'HorizontalAlignment','left')
annotation('textbox', [0.3, 1, 0.4,0], 'string', p.Results.titleText,'HorizontalAlignment','center','Fontsize',20,'interpreter','none');

%Save figure if desired
if ~isempty(p.Results.saveFigure)    
    set(f,'Units','Inches','renderer','painters');
    pos = get(f,'Position');
    set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    
        [fpath,fname,fext] = fileparts(p.Results.saveFigure);
        switch(fext)
            case '.pdf'
                print(f,p.Results.saveFigure,'-dpdf','-r0');
            case '.png'
               	print(f,p.Results.saveFigure,'-dpng','-r0');

        end
end

end

function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)

% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering


if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .1; end
if nargin<5; marg_w = .1; end

if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end

axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;

py = 1-marg_h(2)-axh;

% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        %         ha(ii) = axes('Units','normalized', ...
        %             'Position',[px py axw axh], ...
        %             'XTickLabel','', ...
        %             'YTickLabel','');
        
        %PZH changed to always have tick labels
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh]);
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end