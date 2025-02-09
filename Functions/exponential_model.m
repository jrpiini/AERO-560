function rho = exponential_model(alt)

if (alt <= 1000) && (alt > 900)
    H = 181.05; % km
    rho0 = 5.245e-15; % kg/m^3
    h0 = 900; % km
elseif (alt <= 900) && (alt > 800)
    H = 124.64; % km
    rho0 = 1.170e-14; % kg/m^3
    h0 = 800; % km
elseif (alt <= 800) && (alt > 700)
    H = 88.667; % km
    rho0 = 3.614e-14; % kg/m^3
    h0 = 700; % km
elseif (alt <= 700) && (alt > 600)
    H = 71.835; % km
    rho0 = 1.454e-13; % kg/m^3
    h0 = 600; % km
elseif (alt <= 600) && (alt > 500)
    H = 63.822; % km
    rho0 = 6.967e-13; % kg/m^3
    h0 = 500; % km
elseif (alt <= 500) && (alt > 450)
    H = 60.828; % km
    rho0 = 1.585e-12; % kg/m^3
    h0 = 450; % km
elseif (alt <= 450) && (alt > 400)
    H = 58.515; % km
    rho0 = 3.725e-12; % kg/m^3
    h0 = 400; % km
elseif (alt <= 400) && (alt > 350)
    H = 53.298; % km
    rho0 = 9.518e-12; % kg/m^3
    h0 = 350; % km    
elseif (alt <= 350) && (alt > 300)
    H = 53.628; % km
    rho0 = 2.418e-11; % kg/m^3
    h0 = 300; % km
elseif (alt <= 300) && (alt > 250)
    H = 45.546; % km
    rho0 = 7.248e-11; % kg/m^3
    h0 = 250; % km
elseif (alt <= 250) && (alt > 200)
    H = 37.105; % km
    rho0 = 2.789e-10; % kg/m^3
    h0 = 200; % km
elseif (alt <= 200) && (alt > 180)
    H = 29.740; % km
    rho0 = 5.464e-10; % kg/m^3
    h0 = 180; % km
elseif (alt <= 180) && (alt > 150)
    H = 22.523; % km
    rho0 = 2.070e-9; % kg/m^3
    h0 = 150; % km
elseif (alt <= 150) && (alt > 140)
    H = 16.149; % km
    rho0 = 3.845e-9; % kg/m^3
    h0 = 140; % km
elseif (alt <= 140) && (alt > 130)
    H = 12.636; % km
    rho0 = 8.484e-9; % kg/m^3
    h0 = 130; % km
elseif (alt <= 130) && (alt > 120)
    H = 9.473; % km
    rho0 = 2.438e-8; % kg/m^3
    h0 = 120; % km
elseif (alt <= 120) && (alt > 110)
    H = 7.263; % km
    rho0 = 9.661e-8; % kg/m^3
    h0 = 110; % km
elseif (alt <= 110) && (alt > 100)
    H = 5.877; % km
    rho0 = 5.297e-7; % kg/m^3
    h0 = 100; % km
elseif (alt <= 100) && (alt > 90)
    H = 5.382; % km
    rho0 = 3.396e-6; % kg/m^3
    h0 = 90; % km
elseif (alt <= 90) && (alt > 80)
    H = 5.799; % km
    rho0 = 1.905e-5; % kg/m^3
    h0 = 80; % km
elseif (alt <= 80) && (alt > 70)
    H = 6.549; % km
    rho0 = 8.770e-5; % kg/m^3
    h0 = 70; % km
elseif (alt <= 70) && (alt > 60)
    H = 7.714; % km
    rho0 = 3.206e-4; % kg/m^3
    h0 = 60; % km
elseif (alt <= 60) && (alt > 50)
    H = 8.382; % km
    rho0 = 1.057e-3; % kg/m^3
    h0 = 50; % km
elseif (alt <= 50) && (alt > 40)
    H = 7.554; % km
    rho0 = 3.972e-3; % kg/m^3
    h0 = 40; % km
elseif (alt <= 40) && (alt > 30)
    H = 6.682; % km
    rho0 = 1.774e-2; % kg/m^3
    h0 = 30; % km    
elseif (alt <= 30) && (alt > 25)
    H = 6.349; % km
    rho0 = 3.899e-2; % kg/m^3
    h0 = 25; % km
elseif (alt <= 25) && (alt > 0)
    H = 7.249; % km
    rho0 = 1.225; % kg/m^3
    h0 = 0; % km
else
    H = 268;
    rho0 = 3.019e-15;
    h0 = 1000;
end

rho = rho0*exp(-(alt - h0)/H); % kg/m^3
% rho = rho * 1e9; % kg/m^3 --> kg/km^3

end
















