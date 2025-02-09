% Input: quaternion [eps1; eps2; eps3; eta]
% Output: rotation matrix
function C = quat2C(quat)

C = ((2*(quat(4)^2) - 1)*eye(3)) + (2*quat(1:3)*...
    quat(1:3)') - (2*quat(4)*[0 -quat(3) ...
    quat(2); quat(3) 0 -quat(1); -quat(2) ...
    quat(1) 0]);

end