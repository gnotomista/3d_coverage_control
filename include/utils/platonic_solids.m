function [V, F] = platonic_solids(id, Rho, option_display, face_type)
%% platonic_solids : function to compute and display the five platonic solids.
%
% % About / info
%
% Platonic solids verify Euler formula, V - E + F = 2 with here : 
%
% - V the number of vertices
% - E the number of edges
% - F the number of faces.  
%
%
% Author and support : nicolas.douillet (at) free.fr, 2020.
%
%
% Syntax
%
% platonic_solids(id);
%
% platonic_solids(id, Rho);
%
% platonic_solids(id, Rho, option_display);
%
% platonic_solids(id, Rho, option_display, face_type);
%
% [V, F] = platonic_solids(id, Rho, option_display, face_type);
%
%
% Description
%
% platonic_solids(id) computes the vertex coordinates of the solid #id, its corresponding face set, and displays it.
%
% platonic_solids(id, Rho) allows to change the radius Rho. Default value is Rho = 1; 
%
% platonic_solids(id, Rho, option_display) allows to *enable / disable the display.
%
% platonic_solids(id, Rho, option_display, face_type) uses either default face type (id+2) edges polygon when set
% to 'default' or triangular face type when set to 'triangle'. Since tetrahedron, octahedron, and icosahedron
% default face type are already triangles, this option influences only square (id=2) and dodecahedron (id=5).
%
% [V, F] = platonic_solids(id, Rho, option_display, face_type) also returns vertex and face sets.
%
%
% See also : PLOT::TETRAHEDRON, PLOT::CUBE, PLOT::OCTAHEDRON, PLOT::ICOSAHEDRON, PLOT::DODECAHEDRON, PLOT::SPHERE
%
%
% Input arguments
%
% - id : positive integer scalar double in the set {1;2;3;4;5}.
%        id = 1 -> tetrahedron (fire)
%        id = 2 -> cube (earth)
%        id = 3 -> octahedron (air)
%        id = 4 -> icosahedron (water)
%        id = 5 -> dodecahedron (ether)
%
% - Rho : real positive scalar double, the radius of the sphere circumscribed to the solid.
%
% - option_display : logical *true/false, to enable/disable the display.
%
% - face_type : character string in the set {'default','triangle'}. Case insensitive.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex set. Size(V) = [vertex_nb,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - F = [i1 i2 i3], positive integer matrix double, the face set. Size(T) = [face_nb,3].
%       [ |  |  |]
%
%
% Examples
%
% platonic_solids(1); % tetrahedron
%
% platonic_solids(2,1,true,'triangle'); % triangulated cube
%
% platonic_solids(3,9); % octahedron living in the sphere centred in [0 0 0] and of radius Rho = 9.
%
% [V,T] = platonic_solids(4); % icosahedron
%
% [V,T] = platonic_solids(5); % dodecahedron
% V = V + repmat([1 2 3],[size(V,1),1]); % translate the centre from [0 0 0] to [1 2 3]
%% Input parsing
if nargin < 4
    
    face_type = 'default';
    
    if nargin < 3
        
        option_display = true;
        
        if nargin < 2
            
            Rho = 1;
            assert(nargin > 0,'Not enough input argument.');
            
        else
            
            assert(isreal(Rho) & Rho > 0,'Rho must be a positive real number.');
            
        end
        
    end
    
end
%% Body
phi = 0.5*(1+sqrt(5));
switch id
    
    case 1 % tetrahedron | "fire"              
        
                V = [2*sqrt(2)/3  0          -1/3;...
                     -sqrt(2)/3   sqrt(6)/3  -1/3;...
                     -sqrt(2)/3  -sqrt(6)/3  -1/3;... 
                     0            0          1];
        
                  F = [1 3 2;...
                       1 2 4;...
                       2 3 4;...
                       3 1 4];
                   
                   color = [1 0 0];
        
    case 2 % cube | "earth"
        
        V = (sqrt(3)/3)*[ 1  1  1;...
                         -1  1  1;...
                         -1 -1  1;...
                          1 -1  1;...
                          1  1 -1;...
                         -1  1 -1;...
                         -1 -1 -1;...
                          1 -1 -1];                                            
                      
                      if strcmpi(face_type,'default')
                          
                          F = [1 2 3 4;...
                               8 7 6 5;...
                               1 4 8 5;...
                               2 1 5 6;...
                               3 2 6 7;...
                               4 3 7 8];
                          
                      elseif strcmpi(face_type,'triangle')
                          
                        F = [1 2 4;...
                             2 3 4;...
                             1 4 5;...
                             1 5 2;...
                             2 7 3;...
                             4 3 7;...
                             5 4 8;...
                             2 5 6;...
                             2 6 7;...
                             4 7 8;...
                             5 8 7;...
                             6 5 7];
                          
                      end
                      
                      color = [0 1 0];
        
    case 3 % octahedron | "air"
        
        V = [ 1  0  0;...
              0  1  0;...
             -1  0  0;
              0 -1  0;...
              0  0  1;...
              0  0 -1];
         
         F = [1 2 5;...
              2 3 5;...
              3 4 5;...
              4 1 5;...
              1 6 2;...
              2 6 3;...
              3 6 4;...
              4 6 1];
          
          color = [0 1 1];
        
    case 4 % icosahedron | "water"                                                          
        
        Mrz = [cos(0.4*pi) -sin(0.4*pi) 0;...
               sin(0.4*pi)  cos(0.4*pi) 0;...
               0            0           1];
        
        centre_angle = 2*asin(1/sqrt(phi*sqrt(5)));
        % a = 2/sqrt(phi*sqrt(5)); % edge length
        
        % Icosahedron vertices coordinates
        % 1st equilateral triangle
        V0 = [0 0 1]';
        V1 = [sin(centre_angle) 0 cos(centre_angle)]';
        V2 = Mrz*V1;
        
        % Lower base triangle with /O symetry
        V3 = -V0;
        V4 = -V1;
        V5 = -V2;
        
        % (12) vertices set coordinates vector
        U0 = Mrz*V2;
        U1 = Mrz^2*V2;
        U2 = Mrz^3*V2;
        U3 = Mrz*V5;
        U4 = Mrz^2*V5;
        U5 = Mrz^3*V5;
        
        V = [V0 V1 V2 U0 U1 U2 V3 V4 V5 U3 U4 U5]';                
            
            F = [1 2 3;... % top part
                 1 3 4;...
                 1 4 5;...
                 1 5 6;...
                 1 6 2;...
                
                 8 7 9;... % bottom part
                 7 10 9;
                 7 11 10;...
                 7 12 11;...
                 7 8 12;...
                
                 2 11 3;... % middle belt
                 3 11 12;...
                 3 12 4;...
                 4 12 8;...
                 4 8 5;...
                 5 8 9;...
                 5 9 6;...
                 6 9 10;...
                 2 6 10;...
                 2 10 11];                    
         
         color = [0 0 1];
        
    case 5 % dodecahedron | "ether"   
        
        [V,F] = platonic_solids(4,Rho,false); % from the icosahedron, as its dual polyhedron
        
        V = cell2mat(cellfun(@(r) mean(V(r,:),1),num2cell(F,2),'UniformOutput', false));                
        V = V ./ repmat(sqrt(sum(V.^2,2)),[1,3]);
        
        if strcmpi(face_type,'default')
        
        F = [1 2 3 4 5;... % top pentagon
             10 9 8 7 6;...% bottom pentagon
            
             11 12 13 2 1;... % superior belt
             13 14 15 3 2;...
             15 16 17 4 3;...
             17 18 19 5 4;...
             19 20 11 1 5;...
            
             6 7 18 17 16;...% inferior belt
             7 8 20 19 18;...
             8 9 12 11 20;...
             9 10 14 13 12;...
             10 6 16 15 14];
         
         elseif strcmpi(face_type,'triangle')
            
            F = [1 2 3;...
                 1 3 4;...
                 1 4 5;...
                 10 9 8;...
                 10 8 7;...
                 10 7 6;...
                 11 12 13;...
                 11 13 2;...
                 11 2 1;...
                 13 14 15;...
                 13 15 3;...
                 13 3 2;...
                 15 16 17;...
                 15 17 4;...
                 15 4 3;...
                 17 18 19;...
                 17 19 5;...
                 17 5 4;...
                 19 20 11;...
                 19 11 1;...
                 19 1 5;...
                 6 7 18;...
                 6 18 17;...
                 6 17 16;...
                 7 8 20;...
                 7 20 19;...
                 7 19 18;...
                 8 9 12;...
                 8 12 11;...
                 8 11 20;...
                 9 10 14;...
                 9 14 13;...
                 9 13 12;...
                 10 6 16;...
                 10 16 15;...
                 10 15 14];
            
        end
         
         color = [1 1 0];
        
    otherwise
        
        error('Solid id must be an positive integer in the range |[1; 5]|.');
        
end
V = Rho*V;
if option_display
    
    plot_mesh(V,F,Rho,'none',color,2,false);
    
    % % another display type; uncomment/comment to enable/disable
    % plot_mesh(V,F,Rho,color,[0 0 0],0.5,false);
        
end
end % platonic_solids
%% plot subfunction
function [] = plot_mesh(V, F, Rho, facecolor, edgecolor, linewidth, disp_sphere)
h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);
patch('Faces',F,'Vertices',V,'FaceVertexCData',[0 1 1],'FaceColor',facecolor,'LineWidth',linewidth,'EdgeColor',edgecolor), hold on;        
if disp_sphere
   
    angle_step = pi/180;
    theta = (0:angle_step:pi)';
    phi = (0:angle_step:2*pi);
    
    Xs = Rho*sin(theta)*cos(phi);
    Ys = Rho*sin(theta)*sin(phi);
    Zs = Rho*repmat(cos(theta),[1 length(phi)]);
    
    plot3(Xs(:,1:30:end),  Ys(:,1:30:end),  Zs(:,1:30:end),  '-.', 'Color', [1 1 1]), hold on;
    plot3(Ys(1:30:end,:)', Xs(1:30:end,:)', Zs(1:30:end,:)', '-.', 'Color', [1 1 1]), hold on;
    
end
alpha(0.5);
xlabel('X'), ylabel('Y'), zlabel('Z');
axis equal, axis tight;
ax = gca;
ax.Clipping = 'off';
set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
view(3);
end % plot_mesh
