function vParticles()

load('MCMCOutput','parameters_all','n_params','lowerbounds','upperbounds')
lb = lowerbounds;
ub = upperbounds;

for pn = 1:n_params
    parameter_names{pn} = num2str(pn);
end



% These two options (and uncommenting lines 45 or 46) allows for lines
% where bestPart or modePart, respectively
%load(strcat('bestPart_vector_POP',T))
%load(strcat('modePart_vector_POP',T))

fname='Dotum';	fsize = 12;	lw = 3;

%% Parameter correlations
	% 2D plots
	parameters_all = parameters_all';
	p=1;	NBINS=10;	NCON=300;	nX=100;
		figure('units','centimeters','position',[5,5,25,20],'Name','Particle clouds');hold on;
		for i= 1:n_params
			for j = 1:n_params
				if i > j
					subplot(n_params,n_params,p)
					[Ncount,BinCentre] = hist3([parameters_all(:,i) parameters_all(:,j)],[NBINS NBINS]);
					[c,h]=contourf(BinCentre{1,2},BinCentre{1,1},Ncount,NCON,'linestyle','none');
					colorbar('off');%colorbar;grid off;
					h1=gca;set(h1,'FontName',fname,'FontSize',fsize);
					if i ~= n_params
						set(gca,'XTickLabel',[])
					end
					if j ~= 1
						set(gca,'YTickLabel',[])
					end
					if i == n_params
						xlabel(parameter_names(j),'FontName',fname,'FontSize',fsize);
					end
					if j ==1 
						ylabel(parameter_names(i),'FontName',fname,'FontSize',fsize)
					end
				elseif i == j
					subplot(n_params,n_params,p)
					title(parameter_names(i))
					DX=linspace(min(parameters_all(:,i)),max(parameters_all(:,i)),length(parameters_all)/10);
					DY = density(parameters_all(:,i),DX);
					plot(DX,DY,'linewidth',lw);hold on;
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					%plot(bestParams(i).*ones(1,nX),linspace(0,max(DY),nX),'.-','linewidth',lw)
					%plot(modeParams(i).*ones(1,nX),linspace(0,max(DY),nX),'.-','linewidth',lw)
					%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					h1=gca;set(h1,'FontName',fname,'FontSize',fsize);grid on;
					if i ==n_params
						xlabel(parameter_names(i),'FontName',fname,'FontSize',fsize);
					end
					if j ==1
						ylabel('Prob. Density','FontName',fname,'FontSize',fsize);
					end
					if i ~= n_params
						set(gca,'XTickLabel',[])
					end
					if j ~= 1
						set(gca,'YTickLabel',[])
					end
					%xlabel(parameter_names(i),'FontName',fname,'FontSize',fsize);ylabel('Prob. Density','FontName',fname,'FontSize',fsize);
					axis tight;
					%xlim([lb(j) ub(j)])
				end
				p = p + 1;
			end
		end
		%legend('Param. distr.','Param. best part.')
		
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions below are called by above function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,xo]=density(x,xout,ss,gaus)
%DENSITY  Density estimator using Gaussian kernel
% Y = DENSITY(X,XOUT,S)
% X is the vector of data values.
% The density estimator is evaluated at XOUT points.
% S is a scale factor for the default kernel bandwidth,
% default S = 1.
% Without output arguments the density is plotted.

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.9 $  $Date: 2012/09/27 11:47:35 $

	if nargin<3
	  ss=1;
	end
	if nargin<4
	  gaus=1;
	end

	if nargin<2 | isempty(xout)
	  xmin=min(x); xmax=max(x); xrange=xmax-xmin;
	  if length(x) > 200
		xout=linspace(xmin-0.08*xrange,xmax+0.08*xrange);
	  else
		xout=linspace(mean(x)-4*std(x),mean(x)+4*std(x));
	  end
	end
	y  = zeros(size(xout));
	n  = length(xout);
	nx = length(x);

	%%% see MASS 2nd ed page 181.
	if iqrange(x)<=0
	  s=1.06*std(x)*nx^(-1/5);
	else
	  s=1.06*min(std(x),iqrange(x)/1.34)*nx^(-1/5);
	end
	%  s=1.144*std(x)*nx^(-1/5);
	if ss>0
	  s=ss*s;
	elseif ss<0
	  s = abs(ss);
	end
	if gaus==1
	  % Gaussian kernel
	  for i=1:n
		y(i) = 1/nx*sum(norpf((xout(i)-x)/s))./s;
%         norpf_in = (xout(i)-x)/s;
%         mu = 0;
%         sigma2 = 1;
%         norpf_out=1./sqrt(2*pi*sigma2).*exp(-0.5*(norpf_in-mu).^2 ./sigma2);
%         y(i) = 1/nx*sum(norpf_out)./s;
      end
	elseif gaus==-1
	  % Gamma kernel (still testing)
	  s = s*0.5;
	  for i=1:n
		y(i) = 1/nx*sum(gammapf(xout(i),x./s+1,s));
	  end
	else
	  % triangular kernel
	  s=s*1.2113;
	  for i=1:n
		y(i) = 1/nx*sum(max(0,1-abs(xout(i)-x)/s))./s;
	  end
	end

	if nargout>1
	  xo=xout;
	end

	if nargout==0
	  plot(xout,y)
	  clear y % no output
	end
end

function y=iqrange(x)
	% Interquantile range of each column of x

	% ML 2000

	% Marko Laine <marko.laine@fmi.fi>
	% $Revision: 1.3 $  $Date: 2012/09/27 11:47:37 $

	[n,m]=size(x);
	if n==1
	  x=x';
	  n = m;
	  m = 1;
	end

	x  = sort(x);
	% n  = length(x);
	i1 = floor((n+1)/4); 
	i3 = floor(3/4*(n+1));
	f1 = (n+1)/4-i1; 
	f3 = 3/4*(n+1)-i3;
	q1 = (1-f1).*x(i1,:)+f1.*x(i1+1,:);
	q3 = (1-f3).*x(i3,:)+f3.*x(i3+1,:);
	y  = q3-q1;
end

function y=norpf(x,mu,sigma2)
	% NORPF(x,mu,sigma2)  Normal (Gaussian) density function

	% Marko Laine <marko.laine@fmi.fi>
	% $Revision: 1.4 $  $Date: 2012/09/27 11:47:38 $

	if nargin < 2, mu=0; end
	if nargin < 3, sigma2=1; end
	y=1./sqrt(2*pi*sigma2).*exp(-0.5*(x-mu).^2 ./sigma2);
end