function plotMesh(X,Y,Z)
    %OUTPUT:
    %Display a contour plot of the solution on the mesh   

    %total number of elements in the mesh
    T=delaunay(X,Y);
    
    Netot = length(T);
    
    triX=zeros(Netot,3); %x-coordinates of the vertices of each element 
    triY=zeros(Netot,3); %y-coordinates of the vertices of each element
    triVal=zeros(Netot,3); %solution at the vertices of each element

    %loop over elements 
    for e=1:Netot
        %nodes of the current triangle
        nodes = T(e,1:3); 
        %x- and y-coordinatess of nodes in current triangle
        triX(e,:) = X(nodes,1);
        triY(e,:) = Y(nodes);
        %nodal values in current triangle
        triVal(e,:) = Z(nodes);

    end


    
    %Fill each triangle with color gradient based on solution values at the
    %vertices
    %fill(triX',triY',triVal');
    fill(triX',triY',triVal','LineStyle','None');
    
    %Define and customise color bar  
    cbar=colorbar;
    colormap jet;
    brighten(0.5);
    %caxis([2.8e-3,3.3e-3]);
end