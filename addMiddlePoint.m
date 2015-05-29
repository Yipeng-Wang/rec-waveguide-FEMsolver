function [newElements, newCoordinates, newBoundaries] = addMiddlePoint (oldElements, oldCoordinates, oldBoundaries)

temp = 1;
noOfElements = size(oldElements,1);                 % number of elements
noOfNodes = size(oldCoordinates,1);                 % number of nodes
noOfBoundaries = size(oldBoundaries,1);             % number of boundaries

nodesOnEdge = zeros(3*noOfNodes,5);                 % storing the three nodes on each edge(first three columns)and the coordinates of the middle point(last two columns)

% assign the X,Y coordinates of each edge to the first two columns of
% nodesOnEdge matrix
for h=1:noOfElements
    element=zeros(3,2);                             % storing three node numbers in each element
    
    element(1,1) = oldElements(h,1);           
    element(1,2) = oldElements(h,2);
    element(2,1) = oldElements(h,2);
    element(2,2) = oldElements(h,3);
    element(3,1) = oldElements(h,3);
    element(3,2) = oldElements(h,1);

    for i=1:3
        flag = 1;
        for j=1:temp-1
            if (element(i,1)==nodesOnEdge(j,1)&&element(i,2)==nodesOnEdge(j,2))||(element(i,1)==nodesOnEdge(j,2)&&element(i,2)==nodesOnEdge(j,1))
                flag = 0;
            end
        end
        if (flag == 1)
            nodesOnEdge(temp,1) = element(i,1); 
            nodesOnEdge(temp,2) = element(i,2);
            temp = temp+1;
        end
    end
end

% calculating the number of non-zero rows of nodesOnEdge
temp = 0;

for i=1:3 * noOfNodes                         
    if nodesOnEdge(i,1) ~= 0
        temp = temp+1;
    end
end

nonZeroRows = temp;

% add the node number of the middle point on each edge to the column 3,
% and the coordinates to the column 4 and 5
for i = 1:nonZeroRows                                 
    nodesOnEdge(i,3) = noOfNodes+i;
    nodesOnEdge(i,4) = 0.5*(oldCoordinates(nodesOnEdge(i,1),1)+oldCoordinates(nodesOnEdge(i,2),1));
    nodesOnEdge(i,5) = 0.5*(oldCoordinates(nodesOnEdge(i,1),2)+oldCoordinates(nodesOnEdge(i,2),2));
end

newCoordinates = zeros(noOfNodes+nonZeroRows,2);
newCoordinates(1:noOfNodes,:) = oldCoordinates;
newCoordinates(noOfNodes+1:noOfNodes+nonZeroRows,:) = nodesOnEdge(1:nonZeroRows,4:5);     %creating a new coordinate matrix including all nodes coordinates

% create a new element matrix including 6 columns since each element has 6
% nodes for 2nd order FEM
edgeInElement = zeros(noOfElements,3);              % the matrix describing the relationship between edges and elements

for h=1:noOfElements
    element=zeros(3,2);
    
    element(1,1) = oldElements(h,1);
    element(1,2) = oldElements(h,2);
    element(2,1) = oldElements(h,2);
    element(2,2) = oldElements(h,3);
    element(3,1) = oldElements(h,3);
    element(3,2) = oldElements(h,1);
    
    for i = 1:3
        for j = 1:nonZeroRows
            if (element(i,1)==nodesOnEdge(j,1)&&element(i,2)==nodesOnEdge(j,2))||(element(i,1)==nodesOnEdge(j,2)&&element(i,2)==nodesOnEdge(j,1))
                edgeInElement(h,i)=nodesOnEdge(j,3);
            end
        end
    end
end

newElements = zeros(noOfElements, 6);
newElements(:,1:3) = oldElements;
newElements(:,4:6) = edgeInElement;                 % new element matrix including all elementments and nodes

newBoundaries = zeros(2*noOfBoundaries, 1);
newBoundaries(1:noOfBoundaries,1) = oldBoundaries(:,1);

for h=1:noOfBoundaries
    for i=1:nonZeroRows
        if (oldBoundaries(h,1)==nodesOnEdge(i,1)&&oldBoundaries(h,2)==nodesOnEdge(i,2))||(oldBoundaries(h,1)==nodesOnEdge(i,2)&&oldBoundaries(h,2)==nodesOnEdge(i,1))
            newBoundaries(noOfBoundaries+h,1) = nodesOnEdge(i,3);
        end
    end
end

newBoundaries = sortrows(newBoundaries);

end