function propNodes = heatInterpMarkersToNodes2D(T,pcp,markerX,markerY,nodeX,nodeY,dx,dy)
%interpMarkersToNodes2D interpolates a single material property, A, from
%markers to their surrounding grid nodes using a bilinear interpolation.
    
    % Initialize temp arrays
    num = zeros(length(nodeY),length(nodeX));
    denom = zeros(length(nodeY),length(nodeX)); 
    nx = length(nodeX);
    ny = length(nodeY);
    % Loop through markers
    for m=1:1:length(T)
        % Ignore markers outside of the model domain
        if markerX(m) > max(nodeX) || markerX(m) < min(nodeX) || markerY(m) > max(nodeY) || markerY(m) < min(nodeY)
            continue
        end
        
        % Get index of nearest node to the upper left
        j = fix( (markerX(m) - nodeX(1)) / dx ) + 1;
        i = fix( (markerY(m) - nodeY(1)) / dy ) + 1;
        
        if (j<1)
            j=1;
        end
        if (j>nx-1)
            j=nx-1;
        end
        if (i<1)
            i=1;
        end
        if (i>ny-1)
            i=ny-1;
        end

        % Calculate and and y distances between marker and this node
        distanceX = markerX(m)-nodeX(j);
        distanceY = markerY(m)-nodeY(i);
        % If this marker is within dx of the node, include it in the
        % interpolation
        if (distanceX <= dx) && (distanceY <=dy)
            % Neighboring nodes (L=left, R=right, m=marker)
            % UL (i,j)   UR (i,j+1)
            % 0----------0
            % |          |
            % |          |  m=marker
            % |     m    |
            % 0----------0
            % BL (i+1,j)  BR (i+1,j+1)
            % Determine weights
            wtUL = (1-distanceX/dx) * (1-distanceY/dy);
            wtUR = (distanceX/dx) * (1-distanceY/dy);
            wtBL = (1-distanceX/dx) * (distanceY/dy);
            wtBR = (distanceX/dx) * (distanceY/dy);
            
            % Add terms for UL node
            num(i,j) = num(i,j) + (T(m) * pcp(m) * wtUL); 
            denom(i,j) = denom(i,j) + pcp(m)*wtUL;
            
            % Add terms for UR node
            num(i,j+1) = num(i,j+1) + (T(m) * pcp(m) * wtUR);
            denom(i,j+1) = denom(i,j+1) + pcp(m)*wtUR;
            
            % Add terms for BL node
            num(i+1,j) = num(i+1,j) + (T(m) * pcp(m) * wtBL);
            denom(i+1,j) = denom(i+1,j) + pcp(m)*wtBL;
            
            % Add terms for BR node
            num(i+1,j+1) = num(i+1,j+1) + (T(m) * pcp(m) * wtBR);
            denom(i+1,j+1) = denom(i+1,j+1) + pcp(m)*wtBR;
        end
    end
    % Compute final interpolated values using sum arrays
    propNodes = num ./ denom;

end


