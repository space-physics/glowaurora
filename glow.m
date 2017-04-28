function glow()
% quick demo calling GLOW model from Matlab.
% https://www.scivision.co/matlab-python-user-module-import/

flux = 1;  % [erg ...]
E0 = 1e3; % [eV]

glat = 65.1;
glon = -147.5;
t = '2015-12-13T10'; 

G = py.glowaurora.runglowaurora(flux, E0, t, glat, glon);

ver = xarray2mat(G(1));  
z_km = xarrayind2vector(G(1),'z_km');
lambda = xarrayind2vector(G(1),'wavelength_nm')/10;
sza = G(6); sza = char(sza{1});

plotglow(ver,z_km,lambda,sza,t,glat,glon)
end

function plotglow(ver,z_km,lambda,sza,t,glat,glon)
  figure(1), clf(1)
  ax = axes('nextplot','add');
   
  for i = 1:size(ver,2)
    semilogx(ax,ver(:,i), z_km, 'DisplayName',num2str(lambda(i)))
  end
  
  set(ax,'xscale','log')
  title({[t,' SZA: ',sza,' deg.  (',num2str(glat),',', num2str(glon),')']})
  xlabel('Volume Emission Rate [photons ...]')
  ylabel('altitude [km]')
  ylim([100,400])
  xlim([1e-3,1e3])
  grid('on')
  legend('show')
    
end

function V = xarray2mat(V)
  % convert xarray 2-D array to Matlab matrix

  
V= V{1}.values; 
S = V.shape;
V = cell2mat(cell(V.ravel('F').tolist()));
V = reshape(V,[int64(S{1}), int64(S{2})]);
    
end

function I = xarrayind2vector(V,key)
    
I = cell2mat(cell(V{1}.indexes{key}.values.tolist)); 
    
end