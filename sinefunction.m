function [ z ] = sinefunction( x, amp, freq, phase, offset )
% Sinefunction
% Defines sine function that is described by amplitude, frequency, phase and offset
% Amplitude (the height)
% Frequency (the number of oscillations between 0 and 2*PI)
% Position of the sine wave
% Phase (the starting angle value)
% The constant y-offset, or DC (direct current)

% x = linspace(0.2*pi,10)

z = amp * sin(freq * x + phase) + offset

end