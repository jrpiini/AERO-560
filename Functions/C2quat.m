% Input: rotation matrix
% Output: quaternion [eps1; eps2; eps3; eta]
function quat = C2quat(C)

eta = sqrt(trace(C) + 1)/2;

if eta == 0
    abs_eps_1 = sqrt(0.5*(C(1,1)+1));
    abs_eps_2 = sqrt(0.5*(C(2,2)+1));
    abs_eps_3 = sqrt(0.5*(C(3,3)+1));

    if abs_eps_1 > 0
        eps_1 = abs(abs_eps_1);
        eps_2 = sign(C(1,2)) * abs_eps_2;
        eps_3 = sign(C(1,3)) * abs_eps_3;
    elseif abs_eps_2 > 0
        eps_1 = sign(C(1,2)) * abs_eps_1;
        eps_2 = abs_eps_2;
        eps_3 = sign(C(2,3)) * abs_eps_3;
    else 
        eps_1 = sign(C(1,3)) * abs_eps_1;
        eps_2 = sign(C(2,3)) * abs_eps_2;
        eps_3 = abs_eps_3;
    end

else
    eps_1 = (C(2,3) - C(3,2))/(4*eta);
    eps_2 = (C(3,1) - C(1,3))/(4*eta);
    eps_3 = (C(1,2) - C(2,1))/(4*eta);
end

quat = [eps_1; eps_2; eps_3; eta];

end

