function [ m, t ] = stempread( filename )
%STEMPREAD Read measurement datafile made with STemp sensors.
%   [M, T] = STEMPREAD('FILENAME') reads a measurement file FILENAME.  
%   The result is returned in M, a two-dimensional matrix, where each
%   row is a single measurement point of time. The measurement time points
%   (in seconds) are returned in vector T.

s = csvread(filename);
m = s(:, 3:end);
t = s(:, 1);
end

