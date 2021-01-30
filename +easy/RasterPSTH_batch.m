function RasterPSTH_batch(ephysStruct, clusterIDs, saveDir, events, splitBy, splitByColours, sortBy, window, psthBinWidth, psthSmoothWidth, titleText)
%   easy.RasterPSTH_batch(ephysStruct, clusterIDs, saveDir, events, splitBy, splitByColours, sortBy, window, psthBinWidth, psthSmoothWidth, titleText)
%   Generates Raster/PSTHs for all clusters within the ephys struct, using
%   the function easy.RasterPSTH, and saves the figures into a directory
%   of choice. To use this function you have to specify inputs:
%   ephysStruct - output of loadKSDir
%   clusterIDs - vector of clusterIDs used for plotting. If empty, plot for all
%       clusters contained within the ephysStruct.
%   saveDir - Directory to save figures
%   events - Sx2 cell array of event times (see easy.RasterPSTH documentation)
%   splitBy - Sx2 cell array of splitting conditions (see easy.RasterPSTH documentation)
%   splitByColours - Sx1 cell array, each element contains Dx3 RGB values
%   for each splitting condition. Can use [] if you don't need this.
%   sortBy - Sx2 cell array of sorting conditions (see easy.RasterPSTH documentation)
%   window - [start end] window relative to event times to compute. e.g.
%       [-0.5 +0.5]
%   psthBinWidth - scalar bin width (seconds) used to calculate PSTH
%   psthSmoothWidth - scalar width of window (seconds) used for causal
%       gaussian smoothing
%   titleText - string containing title at the top of each plot. Often
%       useful to use expRef for this.

if isempty(clusterIDs)
    clusterIDs = ephysStruct.cids;
end

for c = 1:length(clusterIDs)
    cluID = clusterIDs(c); %cluID is 0-indexed. Therefore cluID+1 converts to 1-indexing.
    spikeTimes = ephysStruct.st(ephysStruct.clu==cluID);
    
    clusterTitle = {titleText,sprintf('clu%d',cluID)};

    fig_path = fullfile(saveDir, [strjoin(clusterTitle,'_'), '.png']   );

    %get kilosort template used for this cluster
    if cluID+1 <= size(ephysStruct.temps,1) %if clusterID has an associated kilosort templtae
        kilosortTemplate = permute(ephysStruct.temps(cluID+1,:,:),[3 2 1]);
    else %otherwise don't plot a template
        kilosortTemplate = [];
    end
    
    easy.RasterPSTH(spikeTimes, events,...
        'splitBy',splitBy,...
        'sortBy',sortBy,...
        'window',window,...
        'splitByColours',splitByColours,...
        'psthBinWidth',psthBinWidth,...
        'psthSmoothWidth',psthSmoothWidth,...
        'titleText',clusterTitle,...
        'kilosortTemplate',kilosortTemplate,...
        'saveFigure',fig_path);
    
    close gcf;
end

end