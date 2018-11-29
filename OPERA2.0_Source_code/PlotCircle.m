% PlotCircle - plots a circle on your figure
% 
% Call the functions as shown below
% 
% PlotCircle(Column,Row,Radius,N,Color);
% 
% PlotCircle is a function that takes 5 inputs
% 
% Inputs Types
%  ------------
%  Column - Integer, Float
%  Row    - Integer, Float
%  Radius - Integer, Float
%  N      - Integer
%  Color  - Character String
% 
% If your figure will be treated as a Matrix (e.g. image)
% -------------------------------------------------------
% 1. Column - is the column(in the matrix) of the center of the circle (Integer)
% 2. Row    - is the row(in the matrix) of the center of the circle (Integer)
% 3. Radius - is the radius of the required circle (Integer)
% 4. N      - is the number of points that will be used to plot the circle (Integer)
% 5. Color  - is the color of the circle (Character String)
% 
% If your figure will be treated as a Normal Graph
% ------------------------------------------------
% 1. Column - is the co-ordinates of the horizintal axis of the center of the circle (Integer, Float)
% 2. Row    - is the co-ordinates of the vertical axis of the center of the circle (Integer, Float)
% 3. Radius - is the radius of the required circle (Integer, Float)
% 4. N      - is the number of points that will be used to plot the circle (Integer)
% 5. Color  - is the color of the circle (Character String)
% 
% Notes on N: The more you increase N, the more you will get an accurate circle
%             The standard value for N is 256
% 
% Notes on Color: Color is a character string, so you must write the charachter between two ''
% 
%  'b'     blue          
%  'g'     green         
%  'r'     red           
%  'c'     cyan            
%  'm'     magenta       
%  'y'     yellow       
%  'k'     black         
%  'w'     white
%  
% Author: Karim Mansour
% E-mail: karim_mansour@msn.com

function PlotCircle(Column,Row,Radius,N,Color)

if(N<=1)
    error('N must be greater than 1');
end

% if (Color ~='b') && (Color ~='g') && (Color ~= 'r') && (Color ~='c') && (Color ~='m') && (Color ~='y') && (Color ~='k') && (Color ~='w')
%     error('This is not an available color, Please use help PlotCircle to choose an appropriate color');
% end

hold on
t=(0:N)*2*pi/N;
x=Radius*cos(t)+Column;
y=Radius*sin(t)+Row;
plot(x,y,'color',Color);
axis square;