function [tab,I,J]=hist3d(x,y,z,xedges,yedges)

%  hist3d: calculates 3d weighted histogramms
%
% Instead of Matlab function hist3 that calculates 2d histograms giving
% the number of elements, hist3d calculates the sum of weights given 
% in the third vector z. The algorithm was written to be fast and stable.
%
% SYNTAX:
%	
%	[tab,I,J]=hist3d(x,y,z,xedges,yedges);
%	
%	INPUTS:
%
%  		x: size=(N,1) or (1,N)
%
%  		y: size=(N,1) or (1,N)
%
%  		z: size=(N,1) or (1,N)
%
%  		xedges: size=(L,1) or (1,L), must be with increasing values,
%          		works also when it's non-linearly spaced.
%
%  		yedges: size=(M,1) or (1,M), must be with increasing values
%          		works also when it's non-linearly spaced
%          		works also when M~=L (with a size different from the one of "xedges").
%
%	OUTPUTS:
%
%  		  tab: matrix whose coefficients are the sum of weights z of elements given by the coordinates x,y that fall in each cell	
%			   size(tab)=(length(yedges),length(xedges))
%
%			I: size=(N,1) gives in which row of the table falls the element given by the coordinate x,y 
%
%			J: size=(N,1) gives in which column of the table falls the element given by the coordinate x,y
%         
% EXAMPLE1: small number of points easy to visualize
%
%      x=[1.2 1.3  2.2  -1.5];
%      y=[1.2 1.8 -0.5   2.5];
%      z=[1.5 -2.5 10  20];
%      xedges=[-3:1:3];
%      yedges=[-5:1:5];
%      [tab,I,J]=hist3d(x,y,z,xedges,yedges);
%      C=tab;C(C==0)=NaN; %colormap
%      surface(xedges,yedges,tab,C)
%
% EXAMPLE2: 4*10^6 elements, two gaussians product
% 
%     R_moy=18; R_std=5; T_moy=30*pi/180; T_std=40*pi/180;
% 
%     T_val=[-100:0.05:100]'*pi/180;   R_val=[0:0.05:35]';
%     [RI,TI]=meshgrid(R_val,T_val);
% 
%     R=  1./(R_std*sqrt(pi*2))*exp(  -(RI-R_moy).^2/(2*R_std^2)  );
%     T=  1./(T_std*sqrt(pi*2))*exp(  -(TI-T_moy).^2/(2*T_std^2)  );
%     DI=R.*T;
% 
%     [x,y]=pol2cart(TI(:)-pi/2,RI(:));
%     z=DI(:);
% 
%     xedges=[-30:0.5:30]; yedges=[-30:0.5:10];
%     [tab,I,J]=hist3d([x;-x],[y;y],[z;z],xedges,yedges);
%     C=tab;C(C<4)=NaN;
%     surface(xedges,yedges,tab,C,'linestyle','none')
%     grid on
%
%   El mehdi Abbou-ou-cherif (Phd student 2014-2017, Blaise Pascal University, Clermont-Ferrand, France)
%   Email        : midi1990@hotmail.fr
%   Last updated : 24/06/2016 with MATLAB 7.0
%
%   This algorithm was partially inspired from hist2 in Matlab File Exchange
%
%   See also HISTC, ACCUMARRAY, SUB2IND



%% ERROR CHECKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<5
         error('Insufficient number of inputs');
end

if (length(x)~=length(y) | length(x)~=length(z) |  length(y)~=length(z) )
         error('the vectors x,y,z must have the same size');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xn, J] = histc(x,xedges);  %J gives the column indices of elements
[yn, I] = histc(y,yedges);  %I gives the row indices of elements

%% ignoring out of range values
idnonnul=find(J'.*I' ~= 0); 
J_new=J(idnonnul);
I_new=I(idnonnul);
z_new=z(idnonnul);

%% creating the output matrix "A" and filling can be done by these steps:  
% A1(I(1),J(1))=z(1)
% A2(I(2),J(2))=z(2)
% A3(I(3),J(3))=z(3)
% ....
% --> A= A1+A2+A3+....An
% a faster way is calculating the indices of the elements in the matrix and
% summing the values of the duplicated indices:

tab=zeros(length(yedges),length(xedges));
id=sub2ind(size(tab),I_new',J_new'); %indices of the matrix
[new_id, f, v] = unique(id(:)); %find duplicated indices
val_summed = accumarray(v, z_new); %summing the values of duplicated indices
tab(new_id)=val_summed;

