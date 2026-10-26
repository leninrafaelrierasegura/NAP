beta_v = linspace(0.25,1.9999,200);
for m=1:4
  rc = zeros(length(beta_v),m);
  rb = zeros(length(beta_v),m+1);
  factor = zeros(length(beta_v),1);
  for i = 1:length(beta_v)
    f = chebfun(@(x) x.^(beta_v(i)-1),[10^(-(5+m)/2),1]);
    [p,q] = chebpade(f,m,m+1);
    c = poly(p);
    b = poly(q);
    rb(i,:) = roots(b)';
    rc(i,:) = roots(c)';
    factor(i) = c(1)/b(1);
  end
  values = [beta_v' factor rc rb];
  name1 = sprintf('%s%s%s','m',num2str(m),'info.txt');
  name2 = sprintf('%s%s%s','m',num2str(m),'.bin');
  fid = fopen(name1,'w');
  fprintf(fid,'%d %d\n',size(values,1),size(values,2));
  fclose(fid)
  fid = fopen(name2,'w');
  fwrite(fid,values,'double');
  fclose(fid);
end