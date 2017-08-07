function ycoeff= yNormalization(selectCoeffs)
% Transform the system's dependent variables to keep their values closer to 1

%ycoeff= [1e1,1,1,1e-3,1e-2,1e-1,1e-2,1e-4,1,1e-1];
ycoeff= [1 1 1 1 1 1 1 1 1 1];
ycoeff= ycoeff(selectCoeffs);
