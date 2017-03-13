function J = olmstead( R, n )
  %n = 40 ;
  invh2 = (n+1)^2 ;

  C = 0.1 ; B = 2 ;
  u = ones(n,1) / sqrt(1+R) ;
  v = ones(n,1) / sqrt(2+R) ;

  L = -2*speye(n,n) + spdiags( ones(n,2), [-1,1], n, n ) ;
  L = L*invh2 ;

  J = sparse( 2*n, 2*n ) ;
  J(1:2:2*n-1,1:2:2*n-1) = C * L + R * speye(n,n) - spdiags(3*u.^2, [0], n, n ) ;
  J(1:2:2*n-1,2:2:2*n) = (1-C) * L - spdiags(2*v, 0, n, n ) ;

  J(2:2:2*n,1:2:2*n-1) = (1./B)*speye(n,n) ;
  J(2:2:2*n,2:2:2*n) = (1./B) * speye(n,n) ;
end
