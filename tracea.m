%provides the partial trace (over subsystem A) for a given bipartite 
% state on systems AB.
function pb=tracea(rho,da)
db=size(rho,1)/da;
pb=zeros(db,db);
for i=1:da
    v=zeros(da,1);v(i)=1;
    k=kron(v',eye(db));
    pb=pb+k*rho*k';
end
end
