% This example script demonstrates Raster and PSTH plotting. First we
% simulate spiking data from a neuron with known responses to some event.
% Then we use the simulated data to run the RasterPSTH function. 

% Simulate spiking as an inhomogeneous poisson process, with firing rate changes aligned to an event
baselineFiringRate = 2; %hz. 
eventFiringRate = 10; %hz. Firing rate of neuron after event.
eventFiringRateDuration = 0.5; %sec duration of firing rate step caused by event
dt = 0.001; %sec, time step for simulation. Must be less than refractoryPeriod
refractoryPeriod = 0.002; %sec, period after each spike when another spike cannot occur
duration = 10000; %sec duration of simulation
eventTimes = linspace( 1, duration-2, 300)'; %list of event times

%Define lambda (firing rate) parameter baseline
time = 0:dt:duration;
lambda = ones(1,length(time)) * baselineFiringRate * dt; %changing firing rate over the time period

%For each event, set the lambda to the event firing rate
for evt = 1:length(eventTimes)
    idx = eventTimes(evt) <= time & time <= eventTimes(evt)+eventFiringRateDuration;
    lambda(idx) = eventFiringRate*dt;
end

%Generate spike times from poisson process
spikeTimes = [];
for t = 1:length(time)
    %coinflip with p=lambda(t). Concatenate spike time if time since last
    %spike > refractoryPeriod
     if binornd(1,lambda(t))==1 && (isempty(spikeTimes) || (time(t)-spikeTimes(end)) > refractoryPeriod )
        spikeTimes = [spikeTimes; time(t)];
     end
end

%Plot raster and PSTH using default parameters
events = {'event time',eventTimes};
title = sprintf('Neuron simulated with %0.2fHz baseline firing rate\nand %0.2fHz firing rate to event lasting %0.2fsec', baselineFiringRate, eventFiringRate, eventFiringRateDuration);
easy.RasterPSTH(spikeTimes, events, 'titleText', title);

%Plot with different smoothing parameters
easy.RasterPSTH(spikeTimes, events, 'titleText', title, 'psthSmoothWidth', 10/1000);
