function [IntInt] = AMP2IntInt(AMP,sigmaXY,sigmaZ)
%Converts amplitude to integrated intensity based on the shape of the
%function. This is a 3D Gaussian integral

IntInt=AMP*(2^(3/2))*(pi^(3/2))*(sigmaXY^2)*sigmaZ;

end

