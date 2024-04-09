close all
clear all
clc

format shorte % see it all is for boys, leave some to the imagination is for men

%
% Description: this code will help us wrap our head around Schoenberg-Whitney
% theorem.
%

% distribute data as you wish (try different stuff!!)
% f = @( x ) sin( x ) + randn( size( x ) ) / 8; % try this too!!
f = @( x ) sin( x );

'initialise x and t';
a  = 1;
b  = 5;

% try different options here, I refer to lsq to the case Nx > Nt, interp when Nx = Nt
% question: what happens when Nx < Nt??
x_is = 'lsq_hardcoded_bad';%'double_linspace';
switch x_is
  case 'test_interp'
    Nt = 5;
    t = linspace( a, b, Nt );
    Nx = Nt;
    x = t + sqrt( eps ); x( end ) = x( end ) - 2 * sqrt( eps );

  case 'lsq_random' % don't know if good or bad
    Nt = 5;
    t = linspace( a, b, Nt );
    Nx = 8;
    x = sort( rand( 1,Nx ) * ( b - a ) + a );

  case 'lsq_hardcoded_bad'
    Nt = 5;
    t = linspace( a, b, Nt );
    x = [ 1.7340e+00
          1.7558e+00
          2.1853e+00
          2.4739e+00
          2.8037e+00
          3.1880e+00
          3.7471e+00
          3.9788e+00 ];
    Nx = size( x,1 );

  case 'lsq_hardcoded_good'
    Nt = 5;
    t = linspace( a, b, Nt );
    x = [ 1.7340e+00
          1.7558e+00
          2.1853e+00
          2.4739e+00
          2.8037e+00
          3.1880e+00
          3.7471e+00
          3.9788e+00 + 1]; % I add one but you know you can add just enough to make this point fall in the last interval
    Nx = size( x,1 );

  case 'interp_random' % don't know if good or bad
    Nt = 5;
    t = linspace( a, b, Nt );
    Nx = Nt;
    x = sort( rand( 1,Nx ) * ( b - a ) + a );

  case 'interp_hardcoded_bad'
    Nt = 5;
    t = linspace( a, b, Nt );
    Nx = Nt;
    x = [ 1.7340e+00
          2.1853e+00
          2.8037e+00
          3.1880e+00
          3.9788e+00 ];

  case 'interp_hardcoded_good'
    Nt = 5;
    t = linspace( a, b, Nt );
    Nx = Nt;
    x = [ 1.7340e+00
          2.1853e+00
          2.8037e+00
          3.1880e+00
          3.9788e+00 + 1];

  case 'double_linspace'
    Nt = 8;
    t = linspace( a,b, Nt );
    Nx = 60;
    x = linspace( a,b, Nx );

end
x = x(:); % make sure of it

% standard meshing operations (start familiarise with it!)
mesh.Points           = t(:);
mesh.ConnectivityList = [ ( 1 : Nt-1 )', ( 2 : Nt )' ];
   B                  = diff( t( mesh.ConnectivityList ), [], 2 );
detB                  = B;

% initialise H
H = zeros( Nx,Nt );
% then loop
for k = 1 : size( mesh.ConnectivityList,1 )
  % once again we assume element's perspective!
  % H( :,k:k+1 ) = H( :,k:k+1 ) + eval_hats_in_element_k( x, k, mesh, B ); % simpler but less general
  H( :,mesh.ConnectivityList(k,:) ) = H( :,mesh.ConnectivityList(k,:) ) + eval_hats_in_element_k( x, k, mesh, B );
end

%{

H1(x1) H2(x1) H3(x1) ... Hn(x1)
H1(x2) H2(x2) H3(x2) ... Hn(x2)
H1(x3) H2(x3) H3(x3) ... Hn(x3)
H1(x4) H2(x4) H3(x4) ... Hn(x4)
H1(x5) H2(x5) H3(x5) ... Hn(x5)
H1(x6) H2(x6) H3(x6) ... Hn(x6)
  .      .      .         .
  .      .      .         .
  .      .      .         .
H1(xm) H2(xm) H3(xm) ... Hn(xm)

m = Nx
n = Nt

%}

' find coefficients';
if Nx == Nt
  c = H \ f( x(:) ); % can you explain the meaning of c?
end
if Nx > Nt
  c = ( H' * H ) \ ( H' * f( x(:) ) ); % why this?
end
if Nx < Nt
  error('Nx must not be smaller than Nt'); % why?
end

'evaluate result!';
xx = linspace( a,b ); % just for visualization purposes
out = 0;
for k = 1 : size( mesh.ConnectivityList,1 )
  out = out + eval_hats_in_element_k( xx(:), k, mesh, B ) * c( k:k+1 ); % could you explain why?
end
plot( xx,  out, '-', 'Linewidth', 2, 'MarkerSize', 10 ), hold on
plot(  t,0 * t, 'x', 'Linewidth', 2, 'MarkerSize', 10 ), hold on % what's this strange thing I'm doing?
plot(  x, f(x), 'o', 'Linewidth', 2, 'MarkerSize', 10 ), hold on



function y = eval_hats_in_element_k( x, k, mesh, B )
  % why is this so different from what you did in class?
  % is it better? (yes)
  % why? (everybody gangsta in 1D)
  phi = @( xi )[ 1 - xi, xi ];
  y = phi( B( k ) \ ( x - mesh.Points( mesh.ConnectivityList( k,1 ),: ) ) ); % what is this?
  y = y .* ( y >= 0 ) .* ( y <= 1 ); % why am I doing this?
end





%
