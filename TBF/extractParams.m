function [f,t,frequency,halflife] = extractParams(useEig)
%% OVERVIEW

% This function extract interpretable parameters from an eigenvalue (pair).

%% Extract frequency.

% Extract the angle in the complex plane.
ang = atan2(imag(useEig),real(useEig));

% Assign the functional parameter and the frequency.
f = ang/(2*pi);
frequency = 100*f;

%% Extract half-life.

% Extract the amplitude of the complex number.
magnitude = abs(useEig);

% Assign the functional form and the halflife.
t = -log(magnitude);
halflife = (-log(1/2)/t)/100;

end

