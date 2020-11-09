% l_ising t_ising XXZ_3D Heisenberg_2D Heisenberg_2D_X

function [simul,H_1_tensor,H_2_tensor,opt4,d] = models(model,opts)

    p = inputParser;
    addParameter(p,'g',1.05);
    addParameter(p,'J',1);
    addParameter(p,'delta',0.5);
    parse(p,opts)
    
    %change this
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %model = "l_ising";



    opt4.method = "svd";
    opt4.to_matrix = 1; %keep in cell form
    opt4.min_single_N=-1;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    S_x_2 = [0,1;1,0];
    S_y_2 = [0,-1i;1i,0];
    S_z_2 = [1,0;0,-1];

    S_x_3 = 1/sqrt(2)* [0,1,0;1,0,1;0,1,0;];
    S_y_3 = 1/sqrt(2)*[0, -1i,0;1i,0, -1i;0, 1i, 0;];
    S_z_3 = [1,0,0;0,0,0;0,0,1;];

    opt4.single_threshold = 1e8;
    opt4.double_threshold = 1e12;
    
    switch model
        case "random"
            d=2;
           
            mask = [0,1i,1i,1i;
                    0,0,1i,1i;
                    0,0,0,1i;
                    0,0,0,0;];
            
            r1 = (rand(d^2,d^2)-1)*2;
            r2 = mask.*(rand(d^2,d^2)-1)*2;
            
            H_2 = ( r1+r2) + ( r1+r2)' ;
            %H_2 = r1+r1';
            
            
            S = svds(H_2)
            
            H_2 = H_2/(sum(S));
            
            H_2_tensor = permute(reshape(H_2,[d,d,d,d]), [3,4,1,2]);
            
            
            H_1 = (rand(2,2)-1)*2;
            b = 1i*rand;
            H_1(1,2)=b;
            H_1(2,1)= -b;
            
            
%             H_2_tensor= -J*(ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]})...
%                     +0.5*g*   ncon( {S_x_2,eye(d)}, {[-1,-3],[-2,-4]})...
%                     +0.5*g*   ncon( {eye(d),S_x_2}, {[-1,-3],[-2,-4]}));
            H_1_tensor = H_1  ;
            


            simul.title = sprintf("H=random");
        
        case "t_ising"
            d=2;
            J=p.Results.J;
            g= p.Results.g;

            %H_2_tensor =-J*ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            %H_1_tensor =-J*g*S_x_2 ;
            
            H_2_tensor= -J*(ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]})...
                    +0.5*g*   ncon( {S_x_2,eye(d)}, {[-1,-3],[-2,-4]})...
                    +0.5*g*   ncon( {eye(d),S_x_2}, {[-1,-3],[-2,-4]}));
            H_1_tensor =-0.0*eye(2) ;
            


            simul.title = sprintf("H=ZZ + %.3fX",g);
          case "t_ising_2"
            d=2;
            J=p.Results.J;
            g= p.Results.g;

            H_2_tensor =-J*ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            H_1_tensor =-J*g*S_x_2 ;
            

            opt4.single_threshold = 1e10;
            opt4.double_threshold = 1e13;
            
            %opt4.single_threshold = 1e-11; %1e-8;
            %opt4.double_threshold = 1e-8;

            simul.title = sprintf("H=ZZ + %.3fX",g);    
            
        case "l_ising"
            d=2;
            J=1;
            g= p.Results.g;

            H_2_tensor =-J*ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            H_1_tensor =-J*( S_z_2 + g*S_x_2)  ;

       

            simul.title = sprintf("H=ZZ+Z+%.2f X",g);

        case "XXZ_3D"
            d = 3;
            delta= p.Results.delta;

            H_2_tensor =-( ncon( {S_x_3,S_x_3}, {[-1,-3],[-2,-4]})+...
                ncon( {S_y_3,S_y_3}, {[-1,-3],[-2,-4]})+...
                delta*ncon( {S_x_3,S_x_3}, {[-1,-3],[-2,-4]}) ...
                );
            H_1_tensor =  zeros(d) ;

            simul.title = sprintf("H=XX+YY+%.3fZZ (3D)",delta);

        case "Heisenberg_2D"
            d=2;
            H_2_tensor = -ncon( {S_x_2,S_x_2}, {[-1,-3],[-2,-4]})...
                -ncon( {S_y_2,S_y_2}, {[-1,-3],[-2,-4]})...
                -ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            H_1_tensor = zeros(d);

        

            simul.title = sprintf("XX+YY+ZZ (2D)");

         case "Heisenberg_2D_X"
            d=2;
            g= p.Results.g;
            H_2_tensor =- ncon( {S_x_2,S_x_2}, {[-1,-3],[-2,-4]})...
                -ncon( {S_y_2,S_y_2}, {[-1,-3],[-2,-4]})...
                -ncon( {S_z_2,S_z_2}, {[-1,-3],[-2,-4]});
            H_1_tensor =- g*S_x_2;

            simul.title = sprintf("XX+YY+ZZ+%.3fX (2D)",g);


        otherwise
            error("unknown model")

    end

end