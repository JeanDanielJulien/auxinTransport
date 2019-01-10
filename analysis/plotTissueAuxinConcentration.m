
%% Plot
figure;
center = [ mean(vertices(1:2:end)) mean(vertices(2:2:end)) ];

rectangleThickness = 0.15 ;

%% Create one axis
ax1 = axes;

%% Plot the tissue, cell by cell
hold on;
for indexCell=1:cell_number
    
    
    fill( vertices( 1 + 2*cells_vertices( ...
        1+(vertices_number_per_cell(indexCell):(vertices_number_per_cell(indexCell+1)-1)) ...
        ) ) - center(1) + ...
        period(...
        1+2*(vertices_number_per_cell(indexCell):(vertices_number_per_cell(indexCell+1)-1)) ...
        ) * width , ...
        ...
        vertices( 2 + 2*cells_vertices( ...
        1+(vertices_number_per_cell(indexCell):(vertices_number_per_cell(indexCell+1)-1)) ...
        ) ) - center(2) + ...
        period(...
        2+2*(vertices_number_per_cell(indexCell):(vertices_number_per_cell(indexCell+1)-1)) ...
        ) * height , ...
        auxin(indexCell), ...
        'LineWidth',1);
    
end


%% Set the aspect ratio of the figure, to have equal distances along both axis
x_range = max(vertices(1:2:end)) - min(vertices(1:2:end)) ;
y_range = max(vertices(2:2:end)) - min(vertices(2:2:end)) ;
pbaspect([x_range y_range 1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'box','on');
set(gca,'fontsize',20);
%% Second Axis

ax2 = axes;

hold on;
for indexCell=1:cell_number
    
    cellWalls =  (vertices_number_per_cell(indexCell)+1):vertices_number_per_cell(indexCell+1) ;
    cellVertices = cells_vertices( cellWalls ) ;
    cellVerticesPositions = [ ...
        vertices( 1 + 2*cellVertices ) + period(2*(cellWalls-1)+1) * width ; ...
        vertices( 2 + 2*cellVertices ) + period(2*cellWalls) * height ]';
    
    indexPrevious=(vertices_number_per_cell(indexCell+1)-vertices_number_per_cell(indexCell));
    for indexWall=1:(vertices_number_per_cell(indexCell+1)-vertices_number_per_cell(indexCell))
        
        directionOfWall = cellVerticesPositions(indexWall,:) - cellVerticesPositions(indexPrevious,:);
        directionOfWall = directionOfWall / norm(directionOfWall);
        
        directionToBarycenters = [ directionOfWall(2) ,-directionOfWall(1) ];
        
        rectangleToPlot = [ ...
            cellVerticesPositions(indexPrevious,:) - center ; ...
            cellVerticesPositions(indexWall,:) - center ; ...
            cellVerticesPositions(indexWall,:) - rectangleThickness * directionToBarycenters - rectangleThickness * directionOfWall / 2 - center ; ...
            cellVerticesPositions(indexPrevious,:) - rectangleThickness * directionToBarycenters + rectangleThickness * directionOfWall / 2 - center ];
        
        fill(rectangleToPlot(:,1),rectangleToPlot(:,2),pin1_density(cellWalls(indexWall)));
        
        indexPrevious=indexWall;
    end
    
    
    %              fill( cellVerticesPositions(:,1), cellVerticesPositions(:,2),PIN1Density(indexWall + 6*(indexCell-1)));
    
    
end


%% Link them together
linkaxes([ax1,ax2])
%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,flipud(gray))
colormap(ax2,'winter')
%% Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.157 .11 .0675 .815]);
cb2 = colorbar(ax2,'Position',[.8 .11 .0675 .815]);
x_range = max(vertices(1:2:end)) - min(vertices(1:2:end)) ;
y_range = max(vertices(2:2:end)) - min(vertices(2:2:end)) ;
pbaspect([x_range y_range 1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'box','on');
set(gca,'fontsize',20);