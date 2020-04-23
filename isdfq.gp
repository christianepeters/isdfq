\\ author: Christiane Peters, http://christianepeters.wordpress.com
\\ October 2010
\\ run as
\\gp isdfq.gp
log2(x)=log(x*1.)/log(2.)

q = 31
n = 961;
k = 771;
x = floor(k/2);
w = 48;

log2q=log2(q);
mincost=10000000; bestp=0; bestl=0;
{for(p=1,10,\
  Anum=binomial(x,p); Bnum=binomial(k-x,p);
  for(l=1,floor(log(Anum)/log(q)+p*log(q-1)/log(q))+10,\
    if(q==2,
      ops=0.5*(n-k)^2*(n+k)\
          + ((0.5*k-p+1)+(Anum+Bnum)*(q-1)^p)*l\
          + q/(q-1.)*(w-2*p+1)*2*p*(1+(q-2)/(q-1.))*\
            Anum*Bnum*(q-1)^(2*p)/q^l,\
      ops=(n-k)^2*(n+k)\
          + ((0.5*k-p+1)+(Anum+Bnum)*(q-1)^p)*l\
          + q/(q-1.)*(w-2*p+1)*2*p*(1+(q-2)/(q-1.))*\
            Anum*Bnum*(q-1)^(2*p)/q^l;
    );
    prob=Anum*Bnum*binomial(n-k-l,w-2*p)/binomial(n,w);
    cost=log2(ops)+log2(log2q)-log2(prob);
    \\print(p," ",l," ",cost);
    if(cost<mincost,
      mincost=cost;
      bestp=p; bestl=l;
    );
));}
\p 8
cost=mincost; p=bestp; l=bestl;
print("Given q=",q,", n="n,", k=",k,", w=",w);
print("parameters p=",p, " and l=",l," yield 2^",cost," bit ops");
