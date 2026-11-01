a = 0.001;

m = 8;
mtable = [];
for h = 0.001:0.001:0.999
      [q1,p1] = brasil(h,m);
      c = fliplr(poly(p1));
      b = fliplr(poly(q1));
      [r,p,k] = residue(c,b)
      mtable = [mtable;h,r',p',k'];
end


save(strcat('m_brasil',num2str(m),'t.mat'),'mtable')




x = linspace(0,1,1000)';
f = x.^h;
P = poly(p1);
X = fliplr([ones(1000,1) x x.^2 x.^3 x.^4]);
denom = X*P'; 

Q = poly(q1);
X = fliplr([ones(1000,1) x x.^2 x.^3 x.^4]);
num = X*Q';

figure(1);clf
plot(x,f)
hold on
plot(x,num./denom)
