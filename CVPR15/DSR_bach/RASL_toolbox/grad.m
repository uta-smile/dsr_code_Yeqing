function result = grad(x)

gradXY = TVOP*x;

result = sqrt(sum(gradXY.^2,3));