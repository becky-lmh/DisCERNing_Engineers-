function [smearedX, smearedY] = smear(newXPos,newYPos)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

resolution = 0.00004;
roundingValue = 1/resolution;

smearedX = round(newXPos .* roundingValue)./roundingValue;
smearedY = round(newYPos .* roundingValue)./roundingValue;

end