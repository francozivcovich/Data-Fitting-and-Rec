close all
clear all
clc

format shorte

%
% Description: we implement the weighted spline equations system from Lecture 11.
%
% You're welcome to play around with this code!
%

a = 0;
b = pi;
f   = @( x )       sin( 2 * x );
fp  = @( x )   2 * cos( 2 * x );
fpp = @( x ) - 4 * sin( 2 * x );

n = 5;%7;
t = linspace( a, b, n );
y = f( t );
xx = linspace( a,b, 1000 );

% Approximation in $4, this command is the Matlab equivalent to our
%  cubicsplineequation( t, f( t ), 'complete_end', fpa, fpa, [] );
cs = spline( t, [ fp( t(1) ), y, fp( t(n) ) ]); % the complete spline

plot( t, y, 'bo', xx, fnval( cs, xx ), 'r', xx, f( xx ), 'b', 'LineWidth', 2 ),
axis( [ a, b, min( y ) - .1, max( y ) + .1 ] )
hold on
plot( t, 0 * t, 'o' )
title('Approximation to the Function')
%


figure,
csp = fnder( cs ); % little babies derivative evaluation version of what we did in the other code
y = fp( xx );
plot( xx, fnval( csp, xx ), 'r', xx, y, 'b', 'LineWidth', 2 ),
axis( [ a, b, min( y ) - .1, max( y ) + .1 ] )
hold on
plot( t, 0*t, 'o' )
title('Approximation to the first derivative')

figure,
cspp = fnder( csp ); % little babies derivative evaluation version of what we did in the other code
y = fpp( xx );
plot( xx, fnval( cspp, xx ), 'r', xx, y, 'b', 'LineWidth', 2 ),
axis( [ a, b, min( y ) - .1, max( y ) + .1 ] )
hold on
plot( t, 0*t, 'o' )
title('Approximation to the second derivative')











%
