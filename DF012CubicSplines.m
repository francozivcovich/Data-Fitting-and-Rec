close all
clear all
clc

format shorte

%
% Description: we implement the weighted spline equations system from Lecture 11.
%
% You're welcome to play around with this code!
%
f = @( x ) double( abs( x - .5 ) < .1 );
% f = @( x ) sin( 2 * pi * x );


a = 0;
b = 1;
Nt = 8;
t = linspace( a, b, Nt );

fpa = 2 * pi; % try setting whatever here and what complete end conditions do
fpb = 2 * pi; % try setting whatever here and what complete end conditions do

% Approximation in $4:
% - dim( $4 ) = n + 2
% [
%  I have 4 * (n - 1) parameters to determine and I set:
%    (n - 2) continuity condition
%    (n - 2) 1st derivative continuity
%    (n - 2) 2n derivative continuity
% that amount to 4*n - 4 - 3*n + 6 = n + 2.
% ]
% - I prescribe n interpolation conditons
% - I prescribe 2 end conditons (left and right ofc)
% in how many ways can I set end conditions? Many, we see two: complete end and
% natural end conditions.

% Approximation in $4 with prescribed value of 1st derivatives at the extremes
% Book-keeping:
ppc = cubicsplineequation( t, f( t ), 'complete_end', fpa, fpa, [] );
coefs_c = ppc.coefs; % here I do this pp nonsense just to let you the opportunity to see how would I interface with the mkpp, ppval, etc. Matlab routines

% Approximation in $4 with zero value of 2nd derivatives at the extremes (and tgv technique)
ppn = cubicsplineequation( t, f( t ), 'natural_end', NaN, NaN, 1e14 );
coefs_n = ppn.coefs; % here I do this pp nonsense just to let you the opportunity to see how would I interface with the mkpp, ppval, etc. Matlab routines

% Approximation in V4:
% - dim( V4 ) = 2n
% [
%  I have 4 * (n - 1) parameters to determine and I set:
%    (n - 2) continuity condition
%    (n - 2) 1st derivative continuity
% that amount to 4*n - 4 - 2*n + 4 = 2*n.
% ]
% - I prescribe n interpolation conditons
% - I got to prescribe n other conditions, I choose here to prescribe at each
%   knot the value of the 1st derivative (how? by minimising somehow the
%   bending energy or however stupidly I want to)
% Alternatively, I can prescribe n-2 1st derivative condition and at the extremes
% I prescribe natural end conditions.


% Cubic weighted spline interpolation that set END second derivatives equal
% to 0 (natural end conditions) and actual (appr) first derivatives at knots.
coefs_wn = weightedspline( t, f( t ), [], 'natural_end' );
% ppw = mkpp( t, coefs_wn );

% Cubic weighted spline interpolation that set END second derivatives equal
% to actual (appr) first derivative (natural end conditions) and actual (appr)
% derivatives at knots.
coefs_wc = weightedspline( t, f( t ), [], 'complete_end', fpa, fpb );
% ppc = mkpp( t, coefs_wc );


% just for visualization purposes
xx = linspace( a,b, 1e3 );

% y = spleval( t, coefs_s, xxx )
figure,
plot( xx, spleval( t, coefs_c, xx ), '-.', ...
      xx, spleval( t, coefs_n, xx ), '--', 'LineWidth',5 )
hold on
plot( xx,       f(             xx ), ':k', 'LineWidth',2 )
hold on
plot( xx, spleval( t, coefs_wc, xx ), '-.', ...
      xx, spleval( t, coefs_wn, xx ), '--', 'LineWidth',5 )
legend('complete end continuous second derivative', 'natural end continuous second derivative', 'function to be approximated', 'complete_end weighted spline', 'natural_end weighted spline')



% Here I show a chosen spline and its first two derivatives
% Ideally, we want to observe that in $4 we have continuity of 1st and 2nd ders
% while in V4 this is not granted.
coeff = coefs_wc; % change this to see successive derivatives of the selected splines

figure,
plot( xx, spleval( t, coeff(:,1:2) .* [ 6,2 ]  , xx ),  '.', 'LineWidth', 3   ),
hold on
plot( t,0*t,'o' )
legend('Second derivative of my spline', 'Knots')

figure,
plot( xx, spleval( t, coeff(:,1:3) .* [ 3,2,1 ], xx ),  '.', 'LineWidth', 3   ),
hold on
plot( t,0*t,'o' )
legend('First derivative of my spline', 'Knots')

figure,
plot( xx, spleval( t, coeff, xx ),  '.', 'LineWidth', 3   ),
hold on
plot( t,0*t,'o' )
legend('My spline', 'Knots')




%
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

function coefs = weightedspline( t, y, w, end_condition, fpa, fpb )

  Nt  = length( t );
  h  = diff( t(:) );
  dd = diff( y(:) ) ./ h; % length( dd ) = n - 1

  if ( nargin < 3 ) || isempty( w )
    pow = 3;
    w = 1 ./ ( 1 + dd.^2 ).^pow; % Attention: here we are using dd in place of derivative!!
  else
    w = w(:);
  end

  b = 3 * [  dd( 1 );
             w( 2 : Nt-1 ) .* h( 1 : Nt-2 ) .* dd( 2 : Nt-1 ) ...
           + w( 1 : Nt-2 ) .* h( 2 : Nt-1 ) .* dd( 1 : Nt-2 );
             dd( Nt - 1 ) ];
  %
  dsub = [ w( 1 : Nt-2 ) .* h( 2 : Nt-1 ); 1; 0 ];
  d    = 2 * [ 1; w( 1 : Nt-2 ) .* h( 2 : Nt-1 ) + w( 2 : Nt-1 ) .* h( 1 : Nt-2 ); 1 ];
  dsup = [ 0; 1; w( 2 : Nt-1 ) .* h( 1 : Nt-2 ) ];

  A = spdiags( [ dsub, d, dsup ], -1 : 1, Nt, Nt );
  if ( nargin > 3 ) && strcmp( end_condition, 'complete_end' )
    % this ain't symmetric and there is no hope to recover symmetricity for
    % non-equispaced knots so no need to worry (uncommon to have equispaced
    % stuff in real life scenarios)but it might be the case that you expressly
    % set your knots equispaced because you see the opportunity and you're smart
    % and you know you can recover symmetricity. How can you do? Very simple
    % once you see it once (help yourself with the notes).
    A(  1,: ) = 0; A(  1,1  ) = 1;
    A( Nt,: ) = 0; A( Nt,Nt ) = 1;
    b(   1 ) = fpa;
    b( end ) = fpb;
  end

  m = A \ b; % stupid stupid stupid
  % in real life we would at least factorise A = L*U or, for large large systems
  % get some preconditioners with ilu (help lu / help ilu)

  coefs( :,4 ) = y( 1 : Nt-1 );
  coefs( :,3 ) = m( 1 : Nt-1 );
  coefs( :,2 ) = ( 3 * dd - 2 * m( 1 : Nt-1 ) - m( 2 : Nt ) ) ./ h; % that's a funny one... (hint: compute s''(t_j))
  coefs( :,1 ) = ( m( 2 : Nt ) - 2 * dd + m( 1 : Nt-1 ) ) ./ h.^2;

end

function pp = cubicsplineequation( t, y, end_condition, fpa, fpb, tgv )

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

  % Compute RHS
  T = [ t( 1 ); t; t( Nt ) ]; % see slides
  Y = [ y( 1 ); y; y( Nt ) ];
  rhs = diff( Y ) ./ diff( T ); % these are actually 1st divided differences right?
  rhs(   1 ) = fpa;
  rhs( end ) = fpb;
  rhs = 3 * diff( rhs );
  % could've replaced these last 6 lines with
  % 3 * diff( [ fpa; diff( y ) ./ h; fpb ] )

  if strcmp( end_condition, 'natural_end' )
    G(   1,: ) = 0;
    G( end,: ) = 0;
    rhs(   1 ) = 0;
    rhs( end ) = 0;
    if ( nargin > 5 ) && not( isempty( tgv ) ) && not( tgv == 0 )
      G(   1,  1 ) = tgv; G(   1,    2 ) = G(     2,  1 );
      G( end,end ) = tgv; G( end,end-1 ) = G( end-1,end );
    else
      G(   1,  1 ) = 1;
      G( end,end ) = 1;
    end
  end

  % Solve ( 6 * G ) * c = rhs
  % we now behave as if it was nothing but in reality the system should be
  % solved using some factorisation. If symmetricity is preserved the best thing
  % to do is use cholesky (see help chol and help ichol):
  % R = chol( G ); % this means G = R'*R -> inv( G ) = inv( R )* inv( R' )
  % C = R \ ( R' \ rhs ) / 6
  % for very large matrices use ichol instead and feed R, R' as preconditioners
  C = ( 6 * G ) \ rhs; % remember: fpp = 2 * c
  A = y;
  D = diff( 2 * C ) ./ ( 6 * h );
  B = diff( A ) ./ h - C(1:end-1) .* h - D .* h.^2;

  % C and A are as long as t, D and B are shorter
  coefs( :,1 ) = D;
  coefs( :,2 ) = C(1:end-1);
  coefs( :,3 ) = B;
  coefs( :,4 ) = A(1:end-1);

  pp = mkpp( t, coefs ); % pointless matlab structure
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
