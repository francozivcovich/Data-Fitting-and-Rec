close all
clear all
clc

format shorte

%
% Description: we wish to build a software for B-spline interpolation. To do so,
% we choose a k and assemble the matrix R = ( B_{i,k}( x_j ) )_ij. Then, we
% make so that the interpolation conditions s( x_i ) = y_i are met.
% In this code, theory and practice heavily blend as Schoenberg-Withney theorem
% for B-splines rules.
%
% When k = 2 everything is rather standard; moreover in this case if one takes
% the interpolation points x to coincide with the knots t, then R turns out to
% being the identity matrix (I expect you to tell me why).
%
% When k > 2 then everything is more complicated: we need more interpolation
% conditions to satisfy S-W theorem (why and how many?).
%
% Also, what happens when the interpolation conditions are more than strictly
% necessary?
%
% You are welcome to play around with this code and test your understanding.
%

% function to interpolate and interval of interest
a = 0;
b = 5;
f = @( x ) sin( pi * x / 2.5 );

% determine k
k = 4;

% determine knots and interpolation points
Nt = 6;
Nx = 6 + ( k - 2 );

t = linspace( a,b,Nt );
x = linspace( a,b,Nx );
y = f( x ); % interpolation condition

t_ = [ a - (k-1:-1:1), t, b + (1:1:k-1) ];

% assemble R
for i = 1 : ( Nt + ( k - 2 ) )
  Bf{ i } = Bspline( t_, i, k );
end
for i = 1 : length( x )
  for j = 1 : length( Bf )
    R( i,j ) = Bf{ j }( x( i ) ); % here I could've simply written Bspline(t_,i,k)
  end
end
c = R \ y(:); % of course we know better than this, right???


xx = linspace( a,b, 1e2 ); % just for evaluation purposes
bb = naiveBspleval( c, Bf, xx ); % there are better ways to do this but we don't have time

figure,
plot( xx, bb, ':', 'LineWidth', 2 )
hold on
plot( xx, f( xx ), '-' )
plot( t, f( t ), 'xk' )




function y = naiveBspleval( c, Bf, x )
  y = 0;
  for i = 1 : length( c )
    y = y + c( i ) * Bf{ i }( x );
  end
end


function B = Bspline( t,i,k )
  if ( k == 1 )
    B = @( x ) ( x >= t( i ) ) .* ( x < t( i + 1 ) );
    return
  end
  Bl = Bspline( t, i    , k - 1 );
  Br = Bspline( t, i + 1, k - 1 );

  B = @( x )   ( x          - t( i ) ) / ( t( i + k - 1 ) - t( i ) ) .* Bl( x ) ...
             + ( t( i + k ) - x      ) / ( t( i + k ) - t( i + 1 ) ) .* Br( x );

end


%
