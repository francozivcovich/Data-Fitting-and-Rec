close all
clear all
clc

format shorte 

%
% Description: we aim to implement a tool for computing the integrals
%
%  b_r = int_[a,b] f (x)H_r(x) dx and G_rc = int_[a,b] H_r(x)H_c(x) dx
%
% for r,c = 1,2,...,#of knots in [a,b]. We discussed in class how impractical
% and uncommon is to compute these integrals analytically and chosen instead to
% employ some quadrature formula.
%
% In this code
% A. we start by setting up a quadrature routine;
% B. we compute b_r ang G_rc in an element-to-element fashion using quadrature (*);
% C. we refine the assembly to make the code more easily generalisable.
%
% (*) that's called "assembly".
%

%
% Funny things we could be doing next (start thinking about it at home):
%
% 1. produce L_2[f](x) for an f defined on a linear subspace of R^2 (not R)
% 2. produce L_2[f](x) for an f defined on a curve in R^2
% 3. compute the L_2 norm of L_2[f](x)
% 4. compute the L_2 norm of f(x) (you need pull-back/push-forward operators)
% 5. L_2[f](x) interpolation in 2D/3D!
%

'---------------------------------------------------------------------> Part A';

[ x, w ] = gauss1Dquadrature01( 1 );
assert( sum( w ) == 1 )                   % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail
assert( all( ( x <= 1 ) .* ( x >= 0 ) ) ) % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail

[ x, w ] = gauss1Dquadrature01( 2 );
assert( sum( w ) == 1 )                   % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail
assert( all( ( x <= 1 ) .* ( x >= 0 ) ) ) % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail

[ x, w ] = gauss1Dquadrature01( 3 );
assert( sum( w ) == 1 )                   % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail
assert( all( ( x <= 1 ) .* ( x >= 0 ) ) ) % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail

[ x, w ] = gauss1Dquadrature01( 4 );
assert( sum( w ) == 1 )                   % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail
assert( all( ( x <= 1 ) .* ( x >= 0 ) ) ) % these sometimes are called sanity checks: quick and stupid tests that shouldn't fail

% a more convincing test now
a = 0;
b = 1;
f = @( x ) sin( pi * x / 4 );
ref = - ( sqrt( 2 ) / 2 - 1 ) * ( 4 / pi );

quadrature_order = 4;
[ x, w ] = gauss1Dquadrature01( quadrature_order );
res = 0;
for q = 1 : length( w )
  res = res + f( x( q ) ) * w( q );
end
tol = 1e-5;
assert( norm( res - ref ) < tol * norm( res ) )

% if you got here without errors you're good to go



'---------------------------------------------------------------------> Part B';

f = @( x ) sin( 2 * pi * x ); % this is the function we want to approximate with L_2[f](x)
Omega = [ 0,1 ]; % my interval domain over which I want to perform the approximation

Nt = 10; % number of knots t
t = linspace( Omega( 1 ), Omega( 2 ), Nt )'; % <-- you can try to make this weird as you like
h = diff( t ); % measure of each interval
Nk = size( h,1 ); % number of intervals (a.k.a. elements)

quadrature_order = 4; % let's not take risks at this stage, you can try make it less
[ xi, w ] = gauss1Dquadrature01( quadrature_order );

% Level 0: the simplest possible way to implement this (check lecture notes!)
b = zeros( Nt,1 ); % allocate (this is supposed to be a dense array!)
for k = 1 : Nk
  A = 0;
  B = 0;
  for q = 1 : length( w )
    A = A + h( k ) * f( h( k ) * xi( q ) + t( k ) ) * ( 1 - xi( q ) ) * w( q );
    B = B + h( k ) * f( h( k ) * xi( q ) + t( k ) ) * (     xi( q ) ) * w( q );
  end
  b( k     ) = b( k     ) + A;
  b( k + 1 ) = b( k + 1 ) + B;
end
G = spalloc( Nt, Nt, 3 * Nt ); % allocate (this is NOT supposed to be a dense array!)
%
% NB: when you allocate an empty sparse array you should declare how many nonzero
% entries there will be. Of course this most times is not possible to know in
% advance so you say a number that is approximatively correct, 3 * Nt in our case,
% so that worst case you slow down a bit when you exceed pre-allocated memory.
%
for k = 1 : Nk
  A = 0;
  B = 0;
  C = 0;
  D = 0;
  for q = 1 : length( w )
    A = A + h( k ) * ( 1 - xi( q ) ) * ( 1 - xi( q ) ) * w( q );
    B = B + h( k ) * ( 1 - xi( q ) ) * (     xi( q ) ) * w( q );
    C = C + h( k ) * (     xi( q ) ) * ( 1 - xi( q ) ) * w( q );
    D = D + h( k ) * (     xi( q ) ) * (     xi( q ) ) * w( q );
  end
  G( k    , k     ) = G( k    , k     ) + A;
  G( k    , k + 1 ) = G( k    , k + 1 ) + B;
  G( k + 1, k     ) = G( k + 1, k     ) + C;
  G( k + 1, k + 1 ) = G( k + 1, k + 1 ) + D;
end
u = G \ b;
%
% NB: in real life forget about backslash or inverting matrices without factoring
% them, we are doing to just for learning.
%
% See also:
% https://it.mathworks.com/help/matlab/math/iterative-methods-for-linear-systems.html
%

x = linspace( Omega( 1 ), Omega( 2 ), 100 );
figure
plot( x, f(x), '- ', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
plot( t,  u  , '-o', 'Linewidth', 2, 'MarkerSize', 10 ),
title('Part B: f(x) vs L_2[f](x)')

do_I_wanna_test_this = true; % set to false if you plan to get crazy with Nt or your computer will freeze
if do_I_wanna_test_this
  G_anal = diag( [ h(1); h(1:end-1) + h(2:end); h(end) ] / 3 ) + diag( h / 6, +1 ) + diag( h / 6, -1 ); % a manina ok
  tol = 1e-8;
  assert( norm( full( G ) - G_anal ) < tol * norm( G_anal ) )
end


'---------------------------------------------------------------------> Part C';

% level 1: slight vectorisation
phi = @( xi )[ 1 - xi; xi ]; % for any xi you give this will be a nice 2x1 array
b = zeros( Nt,1 );
for k = 1 : Nk
  load = 0; % matlab will promote it to a 2x1 array automatically
  for q = 1 : length( w )
    load = load + ( h( k ) * f( h( k ) * xi( q ) + t( k ) ) * w( q ) ) * phi( xi( q ) );
  end
  b( k : k+1 ) = b( k : k+1 ) + load;
end

G = spalloc( Nt, Nt, 3 * Nt );
for k = 1 : Nk
  gram = 0; % matlab will promote it to a 2x2 array automatically
  for q = 1 : length( w )
    gram = gram + h( k ) * w( q ) * ( phi( xi( q ) ) * phi( xi( q ) )' ); % notice phi(xi)*phi(xi)' is a rank-1 symmetric 2x2 matrix ...
  end
  % ... does this mean that gram is to a rank-1 symmetric 2x2 matrix??? Justify your answer or ask me at some point!
  G( k : k+1, k : k+1 ) = G( k : k+1, k : k+1 ) + gram;
end
u = G \ b;

x = linspace( Omega( 1 ), Omega( 2 ), 100 );
figure
plot( x, f(x), '- ', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
plot( t,  u  , '-o', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
title('Part C: f(x) vs L_2[f](x), 1st improvement')


% level 2: refactoring to have a more FEM-ish look
phi = @( xi )[ 1 - xi; xi ];

% NB: mesh structure is VERY important to understand!
mesh.Points           = t; % this belong to R but it could belong to R^N for what we care
mesh.ConnectivityList = [ ( 1 : Nt-1 )', ( 2 : Nt )' ]; % these are positive integers specifying the indexes of the Points forming the elements!
% In our case mesh.ConnectivityList has 2 columns so we are specifying INTERVALS,
% if there were 3 columns we would be specifying TRIANGLES, 4 columns -> TETRAHEDRA,
% of course mesh.Points has to be embedded in the right dimension (you can't make
% triangles out of 1D mesh.Points but nothing stops you from making
% intervals/triangles out of 3D mesh.Points)
   B                  = diff( t( mesh.ConnectivityList ), [], 2 ); % part of the pull-back/push-forward from/to any element (interval) to the reference element (interval)
detB                  = B; % say something about orientation (I forgot!!!)

b = zeros( Nt,1 );
for k = 1 : Nk
  load = 0;
  for q = 1 : length( w )
    load = load + detB( k ) * f( mesh.Points( mesh.ConnectivityList( k,1 ) ) + B( k ) * xi( q ) ) * phi( xi( q ) ) * w( q );
    % easy: we accumulate in load h( k ) * f( eval at push-forward of xi( q ) ) * phi( xi ) * w( q )
  end
  b( mesh.ConnectivityList( k,: ) ) = b( mesh.ConnectivityList( k,: ) ) + load;
  % NB: mesh.ConnectivityList( k,: ) tells us the global numbering of the Points forming our element (interval)
end

G = spalloc( Nt, Nt, 3 * Nt );
for k = 1 : Nk
  gram = 0;
  for q = 1 : length( w )
    gram = gram + detB( k ) * phi( xi( q ) ) * phi( xi( q ) )' * w( q );
  end
  G( mesh.ConnectivityList( k,: ), mesh.ConnectivityList( k,: ) ) = G( mesh.ConnectivityList( k,: ), mesh.ConnectivityList( k,: ) ) + gram;
  % NB: mesh.ConnectivityList( k,: ) tells us the global numbering of the Points forming our element (interval)
end
u = G \ b;

x = linspace( Omega( 1 ), Omega( 2 ), 100 );
figure
plot( x, f(x), '- ', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
plot( t,  u  , '-o', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
title('Part C: f(x) vs L_2[f](x), 2nd improvement')

figure,
spy( G )
title('spy(G) with UNshuffled knots')
% compare this spy plot with the other: do you remember how important is a narrow
% bandwidth for sparse matrices and how much out of our way we usually go for
% keeping it that way under factorisations? (keywords: ilu ichol, call help on
% Matlab to know more about these).


% level 3: shake it to see if something come loose (robustness)
rp = randperm( Nt );                                        % this is unimportant: we are just creating a messy mesh structure
pr( rp ) = 1:Nt;                                            % this is unimportant: we are just creating a messy mesh structure
mesh.Points           = t( rp );                            % this is unimportant: we are just creating a messy mesh structure
mesh.ConnectivityList = [ pr( 1 : Nt-1 )', pr( 2 : Nt )' ]; % this is unimportant: we are just creating a messy mesh structure
   B                  = diff( t( mesh.ConnectivityList ), [], 2 ); % part of the pull-back/push-forward from/to any element (interval) to the reference element (interval)
detB                  = B; % say something about orientation (I forgot!!!)

% NB: the following is exactly equal to level 2, it means to show that even with
% messed up mesh this machinery we set up together still works fine requiring
% ZERO brain-power after the first implementation ;););););););)
b = zeros( Nt,1 );
for k = 1 : Nk
  load = 0;
  for q = 1 : length( w )
    load = load + detB( k ) * f( mesh.Points( mesh.ConnectivityList( k,1 ) ) + B( k ) * xi( q ) ) * phi( xi( q ) ) * w( q );
  end
  b( mesh.ConnectivityList( k,: ) ) = b( mesh.ConnectivityList( k,: ) ) + load;
end
G = spalloc( Nt, Nt, 3 * Nt );
for k = 1 : Nk
  gram = 0;
  for q = 1 : length( w )
    gram = gram + detB( k ) * phi( xi( q ) ) * phi( xi( q ) )' * w( q );
  end
  G( mesh.ConnectivityList( k,: ), mesh.ConnectivityList( k,: ) ) = G( mesh.ConnectivityList( k,: ), mesh.ConnectivityList( k,: ) ) + gram;
end
u = G \ b;

x = linspace( Omega( 1 ), Omega( 2 ), 100 );
figure
plot( x, f(x), '- ', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
[ ~, id ] = sort( mesh.Points );
plot( mesh.Points( id ),  u( id )  , '-o', 'Linewidth', 2, 'MarkerSize', 10 )
title('Part C: f(x) vs L_2[f](x), 3rd improvement')

figure,
spy( G )
title('spy(G) with SHUFFLED knots')
% compare this spy plot with the other: do you remember how important is a
% narrow band width for sparse matrices and how much we are willing to go out of
% our way for keeping it that way under factorisations? Keywords: ilu, ichol -
% call "help ilu" / "help ichol" on Matlab CommandLine to learn more about these.




%
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
