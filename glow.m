function glow()
%% GLOW model from Matlab.
% https://www.scivision.co/matlab-python-user-module-import/
assert(~verLessThan('matlab', '9.5'), 'Matlab >= R2018b required')

params = py.dict(pyargs('flux', 1, 'E0', 1e3, 'glat', 65.1, 'glon', -147.5, 't0', '2015-12-13T10'));

G = py.glowaurora.runglowaurora(params);

ver = xarray2mat(G{'ver'});  
z_km = xarrayind2vector(G,'z_km');
lambda = xarrayind2vector(G,'wavelength_nm')/10;
sza = G.attrs{'sza'};

plotglow(ver,z_km,lambda,sza, params)
end

function plotglow(ver,z_km,lambda, sza, params)
  figure(1), clf(1)
  ax = axes('nextplot','add');
   
  for i = 1:size(ver,2)
    semilogx(ax,ver(:,i), z_km, 'DisplayName',num2str(lambda(i)))
  end
  
  set(ax,'xscale','log')
  title({[char(params{'t0'}),' SZA: ',num2str(double(sza)),' deg.  (',...
         num2str(params{'glat'}),',',...
         num2str(params{'glon'}),')']})
  xlabel('Volume Emission Rate [photons ...]')
  ylabel('altitude [km]')
  ylim([100,400])
  xlim([1e-3,1e3])
  grid('on')
  legend('show')
    
end

function M = xarray2mat(V)
M = double(py.numpy.asfortranarray(V));
end

function I = xarrayind2vector(V,key)
    
I = cell2mat(cell(V.indexes{key}.values.tolist)); 
    
end