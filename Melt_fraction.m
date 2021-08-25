% Function Melt_fraction()
% This function compute melt fraction (xmelt) and latent heat (hlat)
% at given pressure (ppa), temperature (mtk) and rock type (rock)
% Function returns solution for melt fraction (xmelt) 
% and respective latent heat increment (hlat)
function[xmelt,hlat]=Melt_fraction(ppa,mtk,rock)
% 

% Calculate melt fraction using marker type
P=ppa*1e-6; % Pressure, MPa
tl=0; % Liquidus temperature
switch rock

    % 1 = Sticky air/water, no melting
    case 1
    tl=0;

    % 2 = Sediments: latent heat 300 kJ/kg
    case 2
    % Solidus Temperature
    if (P<1200)
        ts=889+17900/(P+54)+20200/(P+54)^2;
    else
        ts=831+0.06*P;
    end
    % Liquidus temperature
    tl=1262+0.09*P;
    
    % Latent heat
    HL=300000;

    % 3 = Basalt: latent heat 380 kJ/kg
    case 3
    % Solidus Temperature
    if (P<1600) 
        ts=973-70400/(P+354)+77800000/(P+354)^2; 
    else
        ts=935+0.0035*P+0.0000062*P^2;
    end
    % Liquidus temperature
    tl=1423+0.105*P;
    % Latent heat
    HL=380000;

    % 4 = Gabbro: latent heat 380 kJ/kg
    case 4
    % Solidus Temperature
    if (P<1600) 
        ts=973-70400/(P+354)+77800000/(P+354)^2; 
    else
        ts=935+0.0035*P+0.0000062*P^2;
    end
    % Liquidus temperature
    tl=1423+0.105*P;
    % Latent heat
    HL=380000;

    % 5 = Lithospheric mantle (dry): latent heat 400 kJ/kg
    case 5
    % Solidus Temperature
    if (P<10000) 
        ts=1394+0.132899*P-0.000005104*P^2; 
    else
        ts=2212+0.030819*(P-10000);
    end
    % Liquidus temperature
    tl=2073+0.114*P;
    % Latent heat
    HL=400000;

    % 6 = Asthenospheric mantle (dry): latent heat 400 kJ/kg
    case 6
    % Solidus Temperature
    if (P<10000) 
        ts=1394+0.132899*P-0.000005104*P^2; 
    else
        ts=2212+0.030819*(P-10000);
    end
    % Liquidus temperature
    tl=2073+0.114*P;
    % Latent heat
    HL=400000;

    % 7 = Hydrated mantle (wet): latent heat 400 kJ/kg
    case 7
    % Solidus Temperature
    if (P<2400) 
        ts=1240+49800/(P+323); 
    else
        ts=1266-0.0118*P+0.0000035*P^2;
    end
    % Liquidus temperature
    tl=2073+0.114*P;
    % Latent heat
    HL=400000;

    % 8 = Upper continental crust: latent heat 300 kJ/kg
    case 8
    % Solidus Temperature
    if (P<1200) 
        ts=889+17900/(P+54)+20200/(P+54)^2; 
    else
        ts=831+0.06*P;
    end
    % Liquidus temperature
    tl=1262+0.09*P;
    % Latent heat
    HL=300000;

    % 9 = Lower continental crust: latent heat 380 kJ/kg
    case 9
    % Solidus Temperature
    if (P<1600) 
        ts=973-70400/(P+354)+77800000/(P+354)^2; 
    else
        ts=935+0.0035*P+0.0000062*P^2;
    end
    % Liquidus temperature
    tl=1423+0.105*P;
    % Latent heat
    HL=380000;
end
% Melt fraction and latent heat calc, check
xmelt=0;
hlat=0;
if (tl>0)
	% Solidus and liquidus must not entersect
    % in the extrapolation region
    if (ts>tl-100) 
        ts=tl-100;
    end
    % Melt fraction
    xmelt=(mtk-ts)/(tl-ts);
    if (xmelt<0) 
        xmelt=0;
    end
    if (xmelt>1)
        xmelt=1;
    end
	% Latent heat calc 
	hlat=HL*xmelt;
end
	


