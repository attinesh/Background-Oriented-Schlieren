%BOS Optics

f = 100;
SensorSizeH =16.6400;
v_px = 6.45e-3;
FOV = 50;
W = 1066;
ZD = 317.5+120;
ZC = 120+482;
g = 1066;
M2 = (g/f - 1)^-1;
M = SensorSizeH/FOV;
vH = v_px/M2;



Ep = (1+ZC/ZD)*v_px*((g - f)/(g*f));