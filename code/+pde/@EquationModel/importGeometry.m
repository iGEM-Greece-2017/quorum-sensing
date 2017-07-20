function gm = importGeometry(self, geofilename)
% importGeometry - Creates a geometric model from STL data read from a file.
%    G = importGeometry(PDEM,FILENAME) Constructs a geometry object from STL 
%    data imported from a file. FILENAME is a string specifying the pathname
%    and filename of the STL file including the file extension ending in 
%    'stl' or 'STL'. If the file is in the current directory or in a directory 
%    on the MATLAB path, the pathname can be omitted. 
%    The geometry object created is assigned to the Geometry property and a 
%    handle to the geometry is returned in G.  
%
%    Example:  Import a 3D STL geometry and plot.
%      numberOfPDEs = 1;
%      thePde = createpde(numberOfPDEs);
%      gm = importGeometry(thePde,'Block.stl')
%      pdegplot(gm, 'FaceLabels','on')
%
%    See also pdegplot, pde.DiscreteGeometry    

% Copyright 2014-2016 The MathWorks, Inc.

        if ~isempty(self.Geometry)                        
           error(message('pde:pdeModel:noAssemSupport')); 
        end
        
        fid = fopen(geofilename, 'r');

        if (fid == -1)
            dirfile = dir(geofilename);
            if ~isempty(dirfile)
                if (length(dirfile) ~= 1)
                    error(message('pde:importGeometry:fileIsDirectory', geofilename));           
                else
                    error(message('pde:importGeometry:fileReadRights', geofilename));
                end
            else  
                error(message('pde:importGeometry:fileDoesNotExist', geofilename));
            end
        else
            % File exists and may be on the MATLAB path.  
            % Get full pathname and filename.
            geofilename = fopen(fid);
            fclose(fid);
        end

        % Verify the file extension
        [~,~,ext] = fileparts(geofilename);
        if ~strcmpi(ext,'.stl')
           error(message('pde:importGeometry:stlFileExtension', geofilename));
        end

        gm = pde.DiscreteGeometry(geofilename);   
        self.Geometry = gm;
        self.IsTwoD = false;
    end