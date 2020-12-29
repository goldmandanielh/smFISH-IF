function [euclidianDistance] = euclidDistance(pos_1, pos_2)
%Given 2 points, returns euclidan distance

xdist=pos_1(1,1)-pos_2(1,1);
ydist=pos_1(1,2)-pos_2(1,2);

euclidianDistance=sqrt(xdist^2+ydist^2);

end

