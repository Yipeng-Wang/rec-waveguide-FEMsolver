function [Ka, Kb] = derivationOfEquation (Elements, Coordinates)

Ka=zeros(size(Coordinates,1));
Kb=zeros(size(Coordinates,1));

for h=1:size(Elements, 1)   
    x = zeros(6,1);
    y = zeros(6,1);

    for i=1:6
        x(i) = Coordinates(Elements(h,i),1);
        y(i) = Coordinates(Elements(h,i),2);
    end

    A=zeros(6,6);
    B=zeros(6,6);
    
    % sampling points and the weights employed in the project 
    L2 = [0.33333333;0.47014206;0.05961587;0.47014206;0.10128651;0.79742699;0.10128651];
    L3 = [0.33333333;0.47014206;0.47014206;0.05961587;0.10128651;0.10128651;0.79742699];
    weight = [0.22500000;0.13239415;0.13239415;0.13239415;0.12593918;0.12593918;0.12593918];

%     L2=[0.8738219710;0.0630890144;0.0630890144;0.5014265096;0.2492867451;0.2492867451;0.6365024991;0.6365024991;0.3103524510;0.3103524510;0.0531450498;0.0531450498];
%     L3=[0.0630890144;0.0630890144;0.8738219710;0.2492867451;0.2492867451;0.5014265096;0.3103524510;0.0531450498;0.6365024991;0.0531450498;0.3103524510;0.6365024991];
%     weight=[0.0508449063;0.0508449063;0.0508449063;0.1167862757;0.1167862757;0.1167862757;0.0828510756;0.0828510756;0.0828510756;0.0828510756;0.0828510756;0.0828510756];

    
    for i=1:6
        for j=1:6
            for k=1:7
                ksi=L2(k);
                eta=L3(k);
            
                N=zeros(6,1);
                N(1) = (1-ksi-eta)*(1-2*ksi-2*eta);
                N(2) = ksi*(2*ksi-1);
                N(3) = eta*(2*eta-1);
                N(4) = 4*(1-ksi-eta)*ksi;
                N(5) = 4*ksi*eta;
                N(6) = 4*(1-ksi-eta)*eta;
            
                X=0;
                Y=0;
                for t=1:6
                    X=X+N(t,1)*x(t,1);
                    Y=Y+N(t,1)*y(t,1);
                end
            
                N_ksi=zeros(6,1);                        %derivative of N of ksi
                N_ksi(1,1) = -3+4*ksi+4*eta;
                N_ksi(2,1) = 4*ksi-1;
                N_ksi(3,1) = 0;
                N_ksi(4,1) = 4-8*ksi-4*eta;
                N_ksi(5,1) = 4*eta;
                N_ksi(6,1) = -4*eta;
            
                N_eta=zeros(6,1);                        %derivative of N of eta
                N_eta(1,1) = -3+4*ksi+4*eta;
                N_eta(2,1) = 0;
                N_eta(3,1) = 4*eta-1;
                N_eta(4,1) = -4*ksi;
                N_eta(5,1) = 4*ksi;
                N_eta(6,1) = 4-8*eta-4*ksi;
            
                X_ksi = 0;
                Y_ksi = 0;
                X_eta = 0;
                Y_eta = 0;
            
                for t=1:6
                    X_ksi = X_ksi+N_ksi(t,1)*x(t,1);
                    Y_ksi = Y_ksi+N_ksi(t,1)*y(t,1);
                    X_eta = X_eta+N_eta(t,1)*x(t,1);
                    Y_eta = Y_eta+N_eta(t,1)*y(t,1);
                end
                
                J = X_ksi*Y_eta-Y_ksi*X_eta;            % calculate the Jacobian matrix
            
                N_x = zeros(6,1);
                N_y = zeros(6,1);
            
                for t=1:6
                    N_x(t,1) = (Y_eta*N_ksi(t,1)-Y_ksi*N_eta(t,1))/J;
                    N_y(t,1) = (-X_eta*N_ksi(t,1)+X_ksi*N_eta(t,1))/J;
                end
            
                %applying Gaussian quadrature to integrate the function
                %over a triangular element
                A(i,j) = A(i,j)+weight(k,1)*(N_x(i,1)*N_x(j,1)+N_y(i,1)*N_y(j,1))*J*0.5;
                B(i,j) = B(i,j)+weight(k,1)*N(i,1)*N(j,1)*J*0.5;
            end
        end
    end

    for i = 1:6
        for j = 1:6
            Ka(Elements(h,i),Elements(h,j)) = A(i,j)+Ka(Elements(h,i),Elements(h,j));
            Kb(Elements(h,i),Elements(h,j)) = B(i,j)+Kb(Elements(h,i),Elements(h,j));
        end
    end
end

end