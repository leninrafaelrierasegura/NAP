% Isolated audit: never executes or overwrites the repository generator outputs.
% Set audit_output_dir to a NEW diagnostic directory before calling this script.
% Optionally set chebfun_path; the audit does not call savepath.
% Defaults: beta = (0.5:0.1:1.9)/2 and orders 1:10. Override with
% audit_beta_values and audit_orders (use 1:8 for the computed root table).
% Example (use an existing Chebfun installation):
% chebfun_path='/path/to/chebfun'; audit_output_dir='/tmp/new-order-audit';
% audit_orders=1:8;
% run('/path/to/NAP/matlab/audit_rational_generation.m');
if exist('chebfun_path','var'), addpath(chebfun_path); end
if isempty(which('chebpade')), error('Add Chebfun to the path before running this audit.'); end
if ~exist('audit_output_dir','var'), error('An isolated output directory is required.'); end
if ~isfolder(audit_output_dir), mkdir(audit_output_dir); end
if ~exist('audit_beta_values','var'), audit_beta_values = (5:19)/20; end
if ~exist('audit_orders','var'), audit_orders = 1:10; end
validateattributes(audit_beta_values,{'numeric'},{'vector','real','finite','>=',0.25,'<',1});
validateattributes(audit_orders,{'numeric'},{'vector','integer','positive','finite'});
fprintf('MATLAB %s\n', version);
fprintf('Chebfun %s\n', which('chebfun'));
summary_rows = {};
generated = struct('m',{},'beta',{},'factor',{},'rc',{},'rb',{},'p_cheb',{},'q_cheb',{},'p_monomial',{},'q_monomial',{});
for m = audit_orders(:).'
  for beta = audit_beta_values(:).'
    fprintf('Auditing m=%d beta=%.2f\n',m,beta);
    lastwarn('');
    elapsed = tic;
    try
      lower = 10^(-(5+m)/2);
      f = chebfun(@(x) x.^(beta-1), [lower,1]);
      [p,q] = chebpade(f,m,m+1);
      c = poly(p); b = poly(q);
      rc = roots(c); rb = roots(b);
      factor = c(1)/b(1);
      max_imag = max(abs(imag([rc;rb])));
      sorted_c = sort(real(rc)); sorted_b = sort(real(rb));
      interlaced = length(rc)==m && length(rb)==m+1 && ...
        all(sorted_b(1:end-1)<sorted_c) && all(sorted_c<sorted_b(2:end));
      negative = all(real([rc;rb])<0);
      real_roots = all(abs(imag([rc;rb]))<=1e-12*max(1,abs(real([rc;rb]))));
      x = logspace(log10(lower),0,2000)';
      rational_direct = p(x)./q(x);
      rational_roots = factor*prod(x-rc.',2)./prod(x-rb.',2);
      conversion_error = max(abs(rational_roots-rational_direct)./max(1,abs(rational_direct)));
      approximation_error = max(abs(rational_direct./(x.^(beta-1))-1));
      root_approximation_error = max(abs(rational_roots./(x.^(beta-1))-1));
      active_x = 16./(16+[(3*pi/2)^2;(2*pi)^2]);
      direct_active = p(active_x)./q(active_x);
      root_active = factor*prod(active_x-rc.',2)./prod(active_x-rb.',2);
      active_rate_error = max(abs(active_x.^(beta-1)./direct_active-1));
      exported_active_rate_error = max(abs(active_x.^(beta-1)./root_active-1));
      % Compare Chebyshev colleague roots with monomial companion roots.
      rc_cheb = roots(p,'all'); rb_cheb = roots(q,'all');
      cheb_conversion = NaN; root_difference = NaN;
      if length(rc_cheb)==m && length(rb_cheb)==m+1
        cheb_reconstruction = factor*prod(active_x-rc_cheb.',2)./prod(active_x-rb_cheb.',2);
        cheb_conversion = max(abs(cheb_reconstruction-direct_active)./max(1,abs(direct_active)));
        root_difference = max([abs(sort(rc)-sort(rc_cheb))./max(1,abs(sort(rc))); ...
                               abs(sort(rb)-sort(rb_cheb))./max(1,abs(sort(rb)))]);
      end
      [warning_message,warning_id] = lastwarn;
      summary_rows(end+1,:) = {m,beta,lower,length(f),length(rc),length(rb),max_imag,...
        real_roots,negative,interlaced,conversion_error,approximation_error,root_approximation_error,...
        active_rate_error,exported_active_rate_error,cheb_conversion,root_difference,...
        toc(elapsed),warning_id,warning_message,''};
      generated(end+1) = struct('m',m,'beta',beta,'factor',factor,'rc',rc,'rb',rb,'p_cheb',chebcoeffs(p),'q_cheb',chebcoeffs(q),'p_monomial',c,'q_monomial',b);
    catch ME
      summary_rows(end+1,:) = {m,beta,NaN,NaN,NaN,NaN,NaN,false,false,false,...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,toc(elapsed),'','',ME.message};
      fprintf('FAILED m=%d beta=%.2f: %s\n',m,beta,ME.message);
    end
  end
end
names = {'m','beta','interval_lower','chebfun_length','numerator_degree','denominator_degree',...
  'maximum_root_imaginary','roots_numerically_real','roots_negative','roots_interlaced',...
  'root_conversion_relative_error','approximation_relative_error_on_generator_interval','root_approximation_relative_error_on_generator_interval',...
  'active_rate_error_direct','active_rate_error_exported','chebyshev_root_active_conversion_error',...
  'monomial_vs_chebyshev_roots_relative_difference','seconds','warning_id','warning_message','error'};
result = cell2table(summary_rows,'VariableNames',names);
writetable(result,fullfile(audit_output_dir,'matlab_generation_audit.csv'));
save(fullfile(audit_output_dir,'isolated_generated_coefficients.mat'),'generated');
% Export an audit-only real CSV with explicit checks; original files are untouched.
coefficient_rows = {};
for k = 1:length(generated)
  g=generated(k);
  if ~all(abs(imag([g.rc;g.rb]))<=1e-12*max(1,abs(real([g.rc;g.rb])))), continue; end
  coefficient_rows(end+1,:)={g.m,g.beta,'factor',0,g.factor};
  for j=1:length(g.rc), coefficient_rows(end+1,:)={g.m,g.beta,'rc',j,real(g.rc(j))}; end
  for j=1:length(g.rb), coefficient_rows(end+1,:)={g.m,g.beta,'rb',j,real(g.rb(j))}; end
end
coefficients=cell2table(coefficient_rows,'VariableNames',{'m','beta','kind','index','value'});
% Format beta with two decimal places; keep coefficient values at round-trip precision.
fid=fopen(fullfile(audit_output_dir,'isolated_generated_coefficients.csv'),'w');
cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'m,beta,kind,index,value\n');
for k=1:height(coefficients)
  fprintf(fid,'%d,%.2f,%s,%d,%.17g\n',coefficients.m(k),coefficients.beta(k),...
    coefficients.kind{k},coefficients.index(k),coefficients.value(k));
end
clear cleaner;
disp(result(:,{'m','beta','roots_numerically_real','roots_negative','roots_interlaced','active_rate_error_direct','active_rate_error_exported','error'}));

% Export polynomial values at round-trip precision and beta with two decimal places.
fid=fopen(fullfile(audit_output_dir,'polynomial_representations.csv'),'w');
cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'m,beta,kind,index,value\n');
for k=1:length(generated)
  g=generated(k);
  for kind={'p_cheb','q_cheb','p_monomial','q_monomial','rc','rb','factor'}
    field=kind{1}; values=g.(field);
    for j=1:length(values)
      if ~isreal(values(j)), continue; end
      fprintf(fid,'%d,%.2f,%s,%d,%.17g\n',g.m,g.beta,field,j,values(j));
    end
  end
end
clear cleaner;
