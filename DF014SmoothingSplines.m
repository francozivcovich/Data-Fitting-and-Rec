close all
clear all
clc

format shorte

%
% Description: we assemble a smoothing spline for some noisy data following a
% certain law!
% Remember:
%      - p = 0: full smoothing (it gona crash, set it to 0.00001 if you wanto);
%      - p = 1: full fitting.
%

p = .5; % <-------------------------------------------------------------------- play with this

a = -5;
b =  5;
Nt = 21; % let this be a bit high so that you can average away the error
t = linspace( a, b, Nt );

% f = @( x ) 10 - x.^2 / 5;
% sigma =  sqrt( 1 : Nt ) / 10;
% noise = randn( 1 , Nt ) .* sigma;
% y = f( t ) + noise;

% funnier option, feel free to play with
f = @( x ) 1 ./ ( 1 + x.^2 );
sigma = ones( 1,Nt ) / 5;
y = f( t ) + 2 * ( rand( size( t ) ) - .5 ) / 10;



% qui iniziamo a fare sul serio
t = t(:);
y = y(:);

h  = diff( t ); % measure of each interval
Nt = size( t,1 ); % number of knots
Nk = size( h,1 ); % number of intervals (a.k.a. elements)

% Grownups' Gram matrix quadrature
mesh.Points           = t;
mesh.ConnectivityList = [ ( 1 : Nt-1 )', ( 2 : Nt )' ];
   B                  = diff( t( mesh.ConnectivityList ), [], 2 );
detB                  = B;
quadrature_order = 4;
[ xi, w ] = gauss1Dquadrature01( quadrature_order );
phi = @( xi )[ 1 - xi; xi ];
G = spalloc( Nt, Nt, 3 * Nt );
for k = 1 : Nk
  gram = 0;
  for q = 1 : length( w )
    gram = gram + detB( k ) * phi( xi( q ) ) * phi( xi( q ) )' * w( q );
  end
  G( mesh.ConnectivityList( k,: ), mesh.ConnectivityList( k,: ) ) = G( mesh.ConnectivityList( k,: ), mesh.ConnectivityList( k,: ) ) + gram;
end

% the more I look at it the more I am pissed that I didn't come up with
% Francesca's idea :C
% It's ok anyway you'll tell me about it during the exam!
Q = ( spdiags( ones( Nt - 1 ,2 ) .* [ -1,1 ], 0:1, Nt - 2, Nt - 1 ) * ...
      spdiags(                          1./h,   0, Nt - 1, Nt - 1 ) * ...
      spdiags( ones( Nt,2 )      .* [ -1,1 ], 0:1, Nt - 1, Nt     ) )';
D = diag( sigma );


am_little_baby = true;
if am_little_baby
  am_little_baby = false; % time to grow up kiddo
end
matfun = @( x ) p * ( D \ ( D \ x(:) ) ) + ( 1 - p ) * ( Q * ( G( 2:end-1,2:end-1 ) \ ( Q' * x(:) ) ) );
tol    = 1e-6;
maxit  = 1e2;
A = pcg( matfun, p * ( D \ ( D \ y(:) ) ), tol, maxit, [], [], y(:) ); % don't have factorisation bcause of Matlab
% ( it pisses me off that one has to pay top schei for Matlab license and then
% it doesn't exists an implementation of a matrix-free incomplete Cholesky ).

% now find C from A!
% C = [ 0; G( 2:end-1,2:end-1 ) \ ( Q' * A / 2 ); 0 ]; % for little babies only
R = ichol( G( 2:end-1,2:end-1 ) );
C = pcg( G( 2:end-1,2:end-1 ), Q' * A / 2, tol, maxit, R, R' );
C = [ 0; C; 0 ];

D = diff( 2 * C ) ./ ( 6 * h );
B = diff( A ) ./ h - C(1:end-1) .* h - D .* h.^2;

% C and A are as long as t, D and B are shorter
coefs( :,1 ) = D;
coefs( :,2 ) = C(1:end-1);
coefs( :,3 ) = B;
coefs( :,4 ) = A(1:end-1);

xx = linspace( a, b, 1e3 );
figure,
plot( xx,                 f( xx ), '-k', 'LineWidth',2 )
hold on
plot(  t,                       y,  'o', 'LineWidth',5 )
hold on
% our smoothing spline
plot( xx, spleval( t, coefs, xx ), '--', 'LineWidth',5 )
% matlab's smoothing spline
[ pp,p ] = csaps( t,y );
plot( xx, ppval( pp,xx ), '--', 'LineWidth',5 )
legend( 'function', 'noisy data', 'our smoothing spline', 'matlab smoothing spline')

warning('our''s and Matlab''s smoothing spline won''t ever coincide cause we dont know the parameters it uses under the hood')


function y = spleval( t, coefs, x )
% Evaluate splines defined by knots t and coefficients coefs at locations x.
% It emulates what done in ppval without the weird pp structure.
  y = zeros( size( x ) );
  x = x(:);
  ord = size( coefs,2 ); % spline order
  for k = 1 : length( t ) - 1
    id = find( ( x >= t( k ) ) .* ( x <= t( k + 1 ) ) );
    y( id ) = ( ( x( id ) - t( k ) ).^( ord-1:-1:0 ) ) * coefs( k,: )';
  end
end

function [ x,w ] = gauss1Dquadrature01( n, a, b )
  if nargin < 2
    a = 0;
    b = 1;
  end
  switch n
    case 1
      x = 0;
      w = 2;
    case 2
      x = sqrt(3) \ [ -1, 1 ];
      w = [ 1 1 ];
    case 3
      x = [ - sqrt( 3/5 ), 0, sqrt( 3/5 ) ];
      w = [ 5 8 5 ] / 9;
    case 4
      x = [ - sqrt( 3/7 + 2/7 * sqrt(6/5) ), - sqrt( 3/7 - 2/7 * sqrt(6/5) ), sqrt( 3/7 - 2/7 * sqrt(6/5) ) sqrt( 3/7 + 2/7 * sqrt(6/5) ) ];
      w = ( 18 + [ -1 1 1 -1 ] * sqrt( 30 ) ) / 36;
  end
  x = ( b - a ) / 2 * x + ( a + b ) / 2;
  w = ( b - a ) / 2 * w;

  x = x(:);
  w = w(:);
end












%
