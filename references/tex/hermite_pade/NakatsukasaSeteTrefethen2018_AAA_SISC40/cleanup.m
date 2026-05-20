function [r,pol,res,zer,z,f,w] = cleanup(r,pol,res,zer,z,f,w,Z,F)
m = length(z); M = length(Z);
ii = find(abs(res)<1e-13);            % find negligible residues
ni = length(ii);
if ni == 0, return, end
fprintf('%d Froissart doublets\n',ni)
for j = 1:ni
  azp = abs(z-pol(ii(j)));
  jj = find(azp == min(azp),1);
  z(jj) = []; f(jj) = [];             % remove nearest support points
end
for j = 1:length(z)
  F(Z==z(j)) = []; Z(Z==z(j)) = [];
end
m = m-length(ii);
SF = spdiags(F,0,M-m,M-m);
Sf = diag(f);
C = 1./bsxfun(@minus,Z,z.');
A = SF*C - C*Sf;
[~,~,V] = svd(A,0); w = V(:,m);       % solve least-squares problem again
r = @(zz) feval(@rhandle,zz,z,f,w);   
[pol,res,zer] = prz(r,z,f,w);         % poles, residues, and zeros
