close all
clear all
clc

format shorte

%
% Description: here I take Bos code for mass extrapolation of gridded 2D Bspline
% interpolants at different grids and slowly get rid of the Matlab's builtin
% routines by substituting them with ours.
%
% This code should've been explained in excruciating detail. We didn't have time
% though. So scroll through ways and try understanding what I did at every step.
%
% The magic here reside in the function meBspline, sapiently written to fit
% highly vectorised queries.
%



Nx = 12;
Ny = 16;

ax = - 1; bx = 1;
ay = - 1; by = 1;

x = linspace( ax, bx, Nx );
y = linspace( ay, by, Ny );

[ X,Y ] = meshgrid( x,y );
Z = sin( pi * ( X.^2 + 2*Y.^2 ) ); % the function to be interpolated

Nevalx = 40; % number of evaluation points
Nevaly = 30; % number of evaluation points
u = linspace( ax, bx, Nevalx );
v = linspace( ay, by, Nevaly );
[ U,V ] = meshgrid( u,v );
W = zeros( size( U ) );

ways = {'bos-massive',
        'bos-simply-looped',
        'bos-fully-looped',
        'fra-fully-looped',
        'fra-simply-looped',
        'fra-massive'};
%
what = ways{6};
disp(['Executing: ', what ]);

k = 4;


tic
switch what
  case 'bos-massive'
    % most compact
    W = spline( y, spline( x,Z,u )', v )';
  case 'bos-simply-looped'
    % second most compact
    for i = 1 : Nevalx
      z = spline( x, Z, u(i) ); % the horizontal interpolants evaluated at x = u(i)
      W( :,i ) = spline( y,z,v );
    end
  case 'bos-fully-looped'
    for iu = 1 : Nevalx
      for iv = 1 : Nevaly
        for i = 1 : Ny
          % fissi y e prendi x ( prendi semplicemente x )
          z( i ) = spline( x, Z( i,: ), u( iu ) );
        end
        % fissi x e prendi y ( prendi semplicemente y )
        W( iv,iu ) = spline( y, z, v( iv ) );
      end
    end
  case 'fra-fully-looped'
    for iu = 1 : Nevalx
      for iv = 1 : Nevaly
        for i = 1 : Ny
          % fissi y e prendi x ( prendi semplicemente x )
          z( i ) = meBspline( x, Z( i,: ), u( iu ), k );
        end
        % fissi x e prendi y ( prendi semplicemente y )
        W( iv,iu ) = meBspline( y, z, v( iv ), k );
      end
    end
  case 'fra-simply-looped'
    z = zeros( Nevalx, Ny );
    for i = 1 : Ny
      % fissi y e prendi x ( prendi semplicemente x )
      z( :,i ) = meBspline( x, Z( i,: ), u, k );
    end
    for iu = 1 : Nevalx
      W( :,iu ) = meBspline( y, z( iu,: ), v, k );
    end
  case 'fra-massive'
    z  = meBspline( x, Z, u, k );
    W = meBspline( y, z, v, k );
end
toc



figure,
surf( U,V,W ); % plot reults
hold on

plot3( X(:), Y(:), Z(:), 'o', 'MarkerFaceColor', 'k' ) % plot initial gridded points
title('Bicubic Spline interpolation')
hold off



function z = meBspline( t, y, x, k )
  % t: knots
  % y: values
  % x: evaluation point
  % k: spline order
  Nt = length( t );
  t_ = [ t(1) - (k-1:-1:1), t, t(end) + (1:1:k-1) ];
  for i = 1 : ( Nt + 1 * ( k - 1 ) - 1 )
    Bf{ i } = Bspline( t_, i, k );
  end
  for i = 1 : length( t )
    for j = 1 : length( Bf )
      R( i,j ) = Bf{ j }( t( i ) );
    end
  end
  c = R \ y';
  Bfm = zeros( size( c,1 ), length( x ) );
  for j = 1 : size( c,1 )
    Bfm( j,: ) = Bf{ j }( x );
  end
  z = Bfm' * c;

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
