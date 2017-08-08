% This program calculates the electromagnetic field in a 2:1 rectangular
% homogenous waveguide with 2nd order finite element method, and displays
% the electric field in TM32 mode.
% Author: Yipeng Wang

clear;
tic;

% import the mesh data produced by the MATLAB mesh generator. It only
% can be used with the first order FEM
TotalEle = dlmread('Node Number.DAT');          % global node number matrix
TotalCoo = dlmread('Coordinates.DAT');          % node coordinates
Boundaries = dlmread('Boundary.DAT');           % boundary condition matrix

TotalEle = TotalEle';
TotalEle = TotalEle(:,1:3);
TotalCoo = TotalCoo';
Boundaries = Boundaries';
Boundaries = Boundaries(:,1:2);

% calculate the second order mesh elements from the first order meshes
[TotalEle, TotalCoo, Boundaries]= addMiddlePoint(TotalEle, TotalCoo, Boundaries); 

% calculate the coefficients Ka and Kb used in the generalized eigenvalue
% equation
[Ka, Kb] = derivationOfEquation(TotalEle, TotalCoo);

noOfElements = size(TotalEle,1);
noOfNodes = size (TotalCoo,1);                  % total number of nodes
noOfBoundaries = size(Boundaries,1);            % total number of nodes on the boundary

% TE modes, solving magnetic field 
[D1,V1] = eig(Ka,Kb);                           
eigenValuesTE = zeros(noOfNodes,1);
eigenVectorsTE = zeros(noOfNodes);

for i=1:noOfNodes
    eigenValuesTE(i) = V1(i,i);
end

eigenValuesTE = sortrows(eigenValuesTE);        % sort the eigenvalues in ascending order

% rearrange the eigenvectors according to the order of eigenvalues
for i = 1:noOfNodes
    for j = 1:noOfNodes
        if V1(i,i) == eigenValuesTE(j)
            eigenVectorsTE(:,j) = D1(:,i);
        end
    end
end

copyOfBoundaries = Boundaries;                  % create an auxiliary column vector of boundary nodes

% impose Dirichlet boundary condition
for i=1:noOfBoundaries;                     
    
    boundaryNode = copyOfBoundaries(i);         % read the node number k of the boundary node
   
    Ka(boundaryNode,:) = [];                    % deleting kth row and column of the global matrices
    Kb(boundaryNode,:) = [];
    Ka(:,boundaryNode) = [];
    Kb(:,boundaryNode) = [];
    copyOfBoundaries = copyOfBoundaries-ones(noOfBoundaries,1);
end

% TM modes, solving electric field 
[D,V] = eig(Ka,Kb);

eigenValuesTM = zeros(noOfNodes-noOfBoundaries,1);
eigenVectorsTM = zeros(noOfNodes-noOfBoundaries);

for i=1:noOfNodes-noOfBoundaries
    eigenValuesTM(i) = V(i,i);
end

eigenValuesTM=sortrows(eigenValuesTM);          % sort the eigenvalues in ascending order

% rearrange the eigenvectors according to the order of eigenvalues
for i=1:noOfNodes-noOfBoundaries
    for j=1:noOfNodes-noOfBoundaries
        if V(i,i) == eigenValuesTM(j)
            eigenVectorsTM(:,j) = D(:,i);
        end
    end
end

eigenVectorsWithBoundaryTM = zeros(noOfNodes,noOfNodes-noOfBoundaries);   

%adding the boundary nodes to the TM eigenvector results
k = 1;
for i = 1:noOfNodes
    flag = 0;
    for j=1:noOfBoundaries
        if i == Boundaries(j)
            flag = 1;
        end
    end
    if flag == 0
        eigenVectorsWithBoundaryTM(i,:)=eigenVectorsTM(k,:);
        k = k+1;
    end
end

T1=zeros(4*noOfElements,3);
for i=1:noOfElements
    T1(4*i-3,1) = TotalEle(i,1);
    T1(4*i-3,2) = TotalEle(i,4);
    T1(4*i-3,3) = TotalEle(i,6);
    T1(4*i-2,1) = TotalEle(i,4);
    T1(4*i-2,2) = TotalEle(i,2);
    T1(4*i-2,3) = TotalEle(i,5);
    T1(4*i-1,1) = TotalEle(i,4);
    T1(4*i-1,2) = TotalEle(i,5);
    T1(4*i-1,3) = TotalEle(i,6);
    T1(4*i,1) = TotalEle(i,6);
    T1(4*i,2) = TotalEle(i,5);
    T1(4*i,3) = TotalEle(i,3);
end

x = TotalCoo(:,1);
y = TotalCoo(:,2);
time = toc; 

% display the electric field of TM32 mode
z = -eigenVectorsWithBoundaryTM(:,7);
figure (1)

trisurf(T1,x,y,z);
title('TM32')
