function [q_conj] = quatconjugate(q)
q_conj = q;
q_conj(1:3) = q_conj(1:3)*-1;
end

