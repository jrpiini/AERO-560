% Input: rotation matrix
% Output: euler angles [phi theta psi] in radians
function angles = euler_from_C(C)

phi = atan2(C(2,3), C(3,3)); % rad
theta = -asin(C(1,3)); % rad
psi = atan2(C(1,2), C(1,1)); % rad

angles = [phi; theta; psi];

end

