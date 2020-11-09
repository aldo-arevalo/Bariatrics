function error = rms(y,ym)
% Root mean square error.
% E = RMS(Y,Ym) returns the root mean squared error between two
% signals. When used to validate models, Y contains the measured
% data and Ym is the model output. When Y and Ym are matrices,
% RMS is computed over the corresponding columns.

% (c) Robert Babuska, 1996

error = sqrt(sum((y - ym).^2)/length(y));
