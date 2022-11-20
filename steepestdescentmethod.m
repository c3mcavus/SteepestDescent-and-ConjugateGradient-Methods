% define the function and variables
syms x y % x -> x1 , y -> x2 for given function
f = 3 + (x - 1.5*y)^2 + (y - 2)^2;
fi = inline(f); % use inline function for give a numerical value back

iter = 20; k = 0; % define maximum number of iteration and step

% define sk and pk for next point -> x(k+1) = x(k) + s(k)*p(k)
syms sk pk;
pk = -gradient(f); pki = inline(pk); % use inline function for give a numerical value back

x0 = [-4.5; -3.5]; % start point 

pk_value = pki(x0(1),x0(2)); % find pk value for start point
sk = vpasolve(diff(fi(x0(1) + sk*pk_value(1),x0(2) + sk*pk_value(2))));

tol = 1e-4; % set the tolerans

x1 = [x0(1) + sk*pk_value(1);x0(2) + sk*pk_value(2)]; % next point

% draw the function 3d using meshgrid and mesh matlab functions - define all points
a = linspace(5,-5); b = linspace(5,-5); [a,b] = meshgrid(a,b); c = 3 + (a - 1.5.*b).^2 + (b - 2).^2;
mesh(a,b,c); hold on;

% draw first point - fi(x0(1),x0(2)) -> value of f(-4.5, -3.5)
plot3(x0(1),x0(2),fi(x0(1),x0(2)),'.','Color', [0, 0, 0],'MarkerSize',15)
% draw line for connecting the points
plot3([x0(1),x1(1)],[x0(2),x1(2)],[fi(x0(1),x0(2)),fi(x1(1),x1(2))],'r','Color', [0, 0, 0],'linewidth',3)

% define a loop for finding a minimum point and use some criterias
while x1(1) - x0(1) > tol || norm(fi(x1(1),x1(2)) - fi(x0(1),x0(2))) > tol || k < iter
    % if the condition is true, continue with next point x0 -> x1 x1 -> x2
    x0 = x1;
    % draw the next point
    plot3(x0(1),x0(2),fi(x0(1),x0(2)),'.','Color', [0, 0, 0],'MarkerSize',15)
    % calculate again sk pk value for next point
    pki_solution = pki(x0(1),x0(2));
    syms sk; sk = vpasolve(diff(fi(x0(1) + sk*pki_solution(1),x0(2) + sk*pki_solution(2))));
    x1 = [x0(1) + sk*pki_solution(1);x0(2) + sk*pki_solution(2)]; % find next point 
    plot3([x0(1),x1(1)],[x0(2),x1(2)],[fi(x0(1),x0(2)),fi(x1(1),x1(2))],'r','Color', [0, 0, 0],'linewidth',3)
    % continue next step
    k = k + 1;
end

% display end point (minimum point)
disp(x1);