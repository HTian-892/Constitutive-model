% Mohr-Coulomb criterion
ecc = 0.2761;  %10-3 the maximum strain of the compaction of the void portion of rock
s3 = 20;  %MPa confining stress
s1c = 230.37;  %MPa peak strength of the axial stress-strain curve 
e1c = 13.74;  %10-3 peak strain of the axial stress-strain curve
E = 21.187;  %GPa apparent elastic modulus for the skeleton portion of rocks
u = 0.227;  % apparent Poisson��s ratio for the skeleton portion of rocks
a = 34.11;  % internal fractal angle
c = 40.28;  % cohesion
X = 2*sin(a)/3^0.5/(3-sin(a));
Df = 1-(1-ecc)*(s1c-u*2*s3)/E/(e1c-ecc);
Ff = E*(e1c-ecc)*(s1c-s3-s1c*sin(a)-s3*sin(a))/(s1c-2*u*s3)/(1-ecc);
Zf = norminv(Df,0,1);
S0 = 1/(1-Df)^2/(2*pi)^0.5*exp(-0.5*(Zf^2))*((s1c-s3)-(s1c+s3)*sin(a))
F02 = Ff-Zf*S0