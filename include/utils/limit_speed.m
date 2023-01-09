function pDot = limit_speed(pDot, V_MAX)
for i = 1 : size(pDot, 2)
    if norm(pDot) >= V_MAX
        pDot = pDot / norm(pDot) * V_MAX;
    end
end
end