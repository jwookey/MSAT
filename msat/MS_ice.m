% MS_ICE - Generate different ice fabrics
%
% // Part of MSAT - The Matlab Seismic Anisotropy Toolkit //
%
% Produce elastic tensor for different ice fabrics. 
%
% [C,rh] = MS_ice(...)
%
% Usage: 
%     [ C,rh ] = MS_ice()                    
%     [ C,rh ] = MS_ice('iso')
%         Return elastic tensor for isotropic ice.
%
%     [ C,rh ] = MS_ice('cluster',phi)
%         Return elastic tensor for an ice 'cluster' fabric, with a cluster
%         opening angle of phi. The cluster is around the x3 axis.
%
%     [ C,rh ] = MS_ice('girdle',xi)
%         Return elastic tensor for an ice 'girdle' fabric, with a cluster
%         thickness angle of xi. The cluster is around the x1 axis.
% 
% Reference: Maurel et al. 2015, Proc. R. Soc. A 471:20140988. 
%            doi:10.1098/rspa.2014.0988 

% Copyright (c) 2011-2020, James Wookey and Andrew Walker
% Copyright (c) 2007-2011, James Wookey
% All rights reserved.
% 
% Redistribution and use in source and binary forms, 
% with or without modification, are permitted provided 
% that the following conditions are met:
% 
%    * Redistributions of source code must retain the 
%      above copyright notice, this list of conditions 
%      and the following disclaimer.
%    * Redistributions in binary form must reproduce 
%      the above copyright notice, this list of conditions 
%      and the following disclaimer in the documentation 
%      and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names 
%      of its contributors may be used to endorse or promote 
%      products derived from this software without specific 
%      prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS 
% AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
% PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
% THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
% USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
% OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [C,rh] = MS_ice(varargin)

   ice_mode = 0  ; % isotropic (default)
          % = 1    % cluster
          % = 2    % girdle

   phi = NaN ;
   xi  = NaN ;

%  ** process the optional arguments
   iarg = 1 ;
   while iarg <= (length(varargin))
      switch lower(varargin{iarg})
         case 'iso'
            ice_mode = 0 ;
            iarg = iarg + 1 ;
         case 'cluster'
            ice_mode = 1 ;
            phi = varargin{iarg+1} ;
            iarg = iarg + 2 ;
         case 'girdle'
            ice_mode = 2 ;
            xi = varargin{iarg+1} ;
            iarg = iarg + 2 ;               
         otherwise 
            error(['Unknown option: ' varargin{iarg}]) ;   
      end   
   end

   rh = 917.0 % same for all 

   if ice_mode==0
      C = MS_iso(3.871810,1.966181,rh) ;
   elseif ice_mode==1
      assert((phi<=90 & phi>=0), ...
         'MS:ICE:BAD_PHI_VALUE','phi should be 0-90 degrees');
      C = ice_cluster(phi) ;
   elseif ice_mode==2
      assert((phi<=90 & phi>=0), ...
         'MS:ICE:BAD_XI_VALUE','xi should be 0-90 degrees');      
      C = ice_girdle(xi) ;
   else
      error('MS:ICE:UNSUPPORTED_MODE','Unsupport fabric mode selected');
   end

   return
end

function [CC] = ice_girdle(xi)
% Reference: Maurel et al. 2015, Proc. R. Soc. A 471:20140988. 
%            doi:10.1098/rspa.2014.0988 
% uses hardwired elastic coefficients

% calculate the fabric using the equations
A = 14.06 ;
C = 15.24 ;
L = 3.06  ;
N = 3.455 ;
F = 5.88  ;

s=sind(xi) ;
CC(1,1)=(1/15.) * (A*(15-10*s^2+3*s^4)+3*C*s^4+2*(2*L+F)*(5*s^2-3*s^4)) ;
CC(2,2)=(1/120.) * (A*(45+10*s^2+9*s^4)+3*C*(15-10*s^2+3*s^4)+2*(2*L+F)*(15+10*s^2-9*s^4)) ;
CC(4,4)=(1/120.) * ((A+C-2*F)*(15-10*s^2+3*s^4)+12*L*(5-s^4)+40*N*s^2) ;
CC(5,5)=(1/30.) * ((A+C-2*F)*(5*s^2-3*s^4)+3*L*(5-5*s^2+4*s^4)+5*N*(3-s^2)) ;
CC(1,2)=(1/30.) * (3*A*(5-s^4)+(C-4*L)*(5*s^2-3*s^4)-10*N*(3-s^2)+F*(15-5*s^2+6*s^4)) ;
CC(1,3)=CC(1,2) ;
CC(2,1)=CC(1,2) ;
CC(2,3)=CC(2,2)-2*CC(4,4) ;
CC(3,1)=CC(1,2) ;
CC(3,2)=CC(2,2)-2*CC(4,4) ;
CC(3,3)=CC(2,2) ;
CC(6,6)=CC(5,5) ;

end

function [CC] = ice_cluster(phi)
% Reference: Maurel et al. 2015, Proc. R. Soc. A 471:20140988. 
%            doi:10.1098/rspa.2014.0988 
% uses hardwired elastic coefficients

% calculate the fabric using the equations
A=14.06 ;
C=15.24 ;
L=3.06 ;
N=3.455 ;
F=5.88 ;

X=1+cosd(phi)+cosd(phi)^2 ;
Y=cosd(phi)^3+cosd(phi)^4 ;

CC = zeros(6,6) ;
CC(1,1) = (1/120.)*(A*(45+19*X+9*Y)+3*C*(15-7*X+3*Y)+2*(2*L+F)*(15+X-9*Y)) ;
CC(3,3) = (1/15.)*(A*(15-7*X+3*Y)+3*C*(X+Y) + 2*(2*L+F)*(2*X-3*Y)) ;
CC(4,4) = (1/30.)*((A+C-2*F)*(2*X-3*Y)+3*L*(5-X+4*Y)+5*N*(3-X)) ;
CC(6,6) = (1/120.)*((A+C-2*F)*(15-7*X+3*Y)+12*L*(5-X-Y)+40*N*X) ;
CC(1,3) = (1/30.)*(3*A*(5-X-Y)+(C-4*L)*(2*X-3*Y)-10*N*(3-X)+F*(15+X+6*Y)) ;
CC(1,2) = CC(1,1) - 2.*CC(6,6) ;
CC(2,1) = CC(1,1) - 2.*CC(6,6) ;
CC(2,2) = CC(1,1) ;
CC(2,3) = CC(1,3) ;
CC(3,1) = CC(1,3) ;
CC(3,2) = CC(1,3) ;
CC(5,5) = CC(4,4) ;

end 



