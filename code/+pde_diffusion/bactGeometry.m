function [gd,name]= bactGeometry(corner, size, i)
% Form the 2D geometry of a bacterium

gd= [3,4, corner(1),corner(1),corner(1)+2*size,corner(1)+2*size,...
          corner(2),corner(2)+size,corner(2)+size,corner(2)]';
name= ['b',num2str(i)];
