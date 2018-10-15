function propMarkers = interpNodesToMarkers2D(A,markerX,markerY,nodeX,nodeY,...
                                                dx,dy,bounds)
%interpNodesToMarkers2D interpolates a material property, A, from grid
%nodes to nearby markers using a bilinear interpolation.
%
% Note: bounds argument is a length=4 array containing the bounding indices
% for nodes to include in the interpolation: bounds=[minX,maxX,minY,maxY].
    propMarkers = zeros(length(markerX),1);
    
    minX = bounds(1);
    maxX = bounds(2);
    minY = bounds(3);
    maxY = bounds(4);
    
    for m=1:1:length(markerX)
        % Get index of nearest node to the upper left
        j = fix( (markerX(m) - nodeX(1)) / dx ) + 1;
        i = fix( (markerY(m) - nodeY(1)) / dy ) + 1;
        % Prevent interpolation from ghost nodes (relies on function input
        % "bounds")
        if (j<minX)
            j=minX;
        end
        if (j>maxX)
            j=maxX;
        end
        if (i<minY)
            i=minY;
        end
        if (i>maxY)
            i=maxY;
        end
        
        % Calculate and and y distances between marker and this node
        % NOTE: removed abs() for these
        distanceX = markerX(m)-nodeX(j);
        distanceY = markerY(m)-nodeY(i);
        % Neighboring nodes (L=left, R=right, m=marker)
        % UL (i,j)   UR (i,j+1)
        % 0----------0
        % |          |
        % |          |  m=marker
        % |     m    |
        % 0----------0
        % BL (i+1,j)  BR (i+1,j+1)

        % Get property values at surrounding nodes
        UL = A(i,j);
        UR = A(i,j+1);
        BL = A(i+1,j);
        BR = A(i+1,j+1);
        % Determine weights
        wtUL = (1-distanceX/dx) * (1-distanceY/dy);
        wtUR = (distanceX/dx) * (1-distanceY/dy);
        wtBL = (1-distanceX/dx) * (distanceY/dy);
        wtBR = (distanceX/dx) * (distanceY/dy); 
        
        % Interpolate from bounding nodes to marker using bilinear
        % interpolation
       propMarkers(m) = wtUL*UL + wtUR*UR + wtBL*BL + wtBR*BR;

    end
end