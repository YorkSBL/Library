function wrappedphi = wrap(phi)
% wrappedphi = wrap(phi)
% Un-unwraps the angle phi (in radians)

  wrappedphi = angle(exp(i*phi));